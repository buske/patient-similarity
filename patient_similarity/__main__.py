#!/usr/bin/env python3

"""
Compute phenotypic similarity between all pairs of patients in patients.json.
A number of different similarity measures are supported, specifiable with
the --score option.
"""

__author__ = 'Orion Buske (buske@cs.toronto.edu)'

import sys
import logging

from collections import Counter

from . import HPO, HPOIC, Diseases, Orphanet, Patient

AOOS = {
    'HP:0003577': 1.5,  # Congenital onset
    'HP:0003623': 3.0,  # Neonatal onset
    'HP:0003593': 4.0,  # Infantile onset
    'HP:0011463': 5.0,  # Childhood onset
    'HP:0003621': 6.0,  # Juvenile onset;
    'HP:0003581': 8.0,  # Adult onset

    'HP:0011460': 1.0,  # Embryonal onset
    'HP:0011461': 2.0,  # Fetal onset
    'HP:0011462': 7.0,  # Young adult onset
    'HP:0003596': 8.0,  # Middle age onset
    'HP:0003584': 9.0,  # Late onset
    }


def similarity_breakdown(hpoic, patient1, patient2):
    p1_terms = set(patient1.hp_terms)
    p2_terms = set(patient2.hp_terms)
    assert p1_terms and p2_terms

    logging.info('Comparing patients: {}, {}'.format(patient1.id, patient2.id))
    logging.info('Patient 1 terms and IC')
    for t in p1_terms:
        logging.info('  {:.6f}: {}'.format(hpoic.get_term_ic(t), t))

    logging.info('Patient 2 terms and IC')
    for t in p2_terms:
        logging.info('  {:.6f}: {}'.format(hpoic.get_term_ic(t), t))

    clusters = []
    while p1_terms and p2_terms:
        p1_ancestors = set().union(*[t.ancestors() for t in p1_terms])
        p2_ancestors = set().union(*[t.ancestors() for t in p2_terms])
        common_ancestors = p1_ancestors & p2_ancestors

        # Find max-ic ancestor
        score, group_root = max([(hpoic.get_term_ic(t), t) for t in common_ancestors])

        # Pop all terms with that as ancestor
        matched1 = set([t for t in p1_terms if group_root in t.ancestors()])
        matched2 = set([t for t in p2_terms if group_root in t.ancestors()])
        p1_terms.difference_update(matched1)
        p2_terms.difference_update(matched2)
        clusters.append((score, group_root, matched1, matched2))

    if p1_terms:
        clusters.append((0, None, p1_terms, []))
    elif p2_terms:
        clusters.append((0, None, [], p2_terms))

    return clusters

def get_aoo_penalty(aoo1, aoo2, neutral_dist=2.0, max_dist=5.0):
    """Get penalty between 1 (worst match) and 0 (no penalty) for age of onsets"""
    if aoo1 and aoo2:
        dist = abs(AOOS[aoo1] - AOOS[aoo2])
        return max(min(dist, max_dist) - neutral_dist, 0) / (max_dist - neutral_dist)
    else:
        # If missing either or both AOOs, no penalty
        return 0

def compare_patients(hpoic, patient1, patient2, scores=None, use_aoo=False):
    logging.debug('Comparing patients: {}, {}'.format(patient1.id, patient2.id))

    p1_terms = patient1.hp_terms
    p2_terms = patient2.hp_terms
    assert p1_terms and p2_terms
    out = {}

    if use_aoo:
        out['aoo'] = get_aoo_penalty(patient1.onset, patient2.onset)

    p1_ancestors = patient1.ancestors()
    p2_ancestors = patient2.ancestors()
    common_ancestors = p1_ancestors & p2_ancestors  # min
    all_ancestors = p1_ancestors | p2_ancestors  # max

    def ancestor_counts(terms):
        counts = Counter()
        for term in terms:
            counts.update(term.ancestors())
        return counts

    def mica(t1, t2):
        return max([(hpoic.get_term_ic(t), t) for t in t1.ancestors() & t2.ancestors()])[1]

    def resnik(t1, t2):
        # Resnik = IC(MICA)
        return max([hpoic.get_term_ic(t) for t in t1.ancestors() & t2.ancestors()])

    def lin(t1, t2, res=None):
        # Lin = 2*resnik / (IC1 + IC2)
        res = resnik(t1, t2) if res is None else res
        return 2 * res / (hpoic.get_term_ic(t1) + hpoic.get_term_ic(t2))

    def jc(t1, t2, res=None):
        # distJ&C = IC1 + IC2 - 2*resnik
        res = resnik(t1, t2) if res is None else res
        return 1 / ((hpoic.get_term_ic(t1) + hpoic.get_term_ic(t2) - 2 * res) + 1)

    def jaccard(t1, t2):
        return len(t1.ancestors() & t2.ancestors()) / len(t1.ancestors() | t2.ancestors())

    # Jaccard = fraction of overlapping nodes
    if not scores or 'ui' in scores:
        out['ui'] = len(common_ancestors) / len(all_ancestors)

    if not scores or 'simgic' in scores:
        out['simgic'] = hpoic.information_content(common_ancestors) / hpoic.information_content(all_ancestors)

    if not scores or 'icca' in scores:
        out['icca'] = hpoic.information_content(common_ancestors)

    if not scores or 'jz' in scores or 'nsimgic' in scores:
        p1_ancestor_counts = ancestor_counts(p1_terms)
        p2_ancestor_counts = ancestor_counts(p2_terms)
        common_ancestor_counts = p1_ancestor_counts & p2_ancestor_counts  # min
        all_ancestor_counts = p1_ancestor_counts | p2_ancestor_counts  # max
        if not scores or 'jz' in scores:
            out['simgic_jz'] = 2 * hpoic.counter_ls_information_content(common_ancestor_counts) / (hpoic.information_content(p1_terms) + hpoic.information_content(p2_terms))

        if not scores or 'nsimgic' in scores:
            out['nsimgic'] = hpoic.counter_information_content(common_ancestor_counts) / hpoic.counter_information_content(all_ancestor_counts)


    if not scores or 'jaccard' in scores:
        jaccard_rows = [[jaccard(t1, t2) for t2 in p2_terms] for t1 in p1_terms]
        out['jaccard_best_avg'] = sum([max(row) for row in jaccard_rows]) / len(jaccard_rows)
        out['jaccard_avg'] = sum([sum(row) for row in jaccard_rows]) / (len(jaccard_rows) * len(jaccard_rows[0]))

    if not scores or ('resnik' in scores or 'owlsim' in scores):
        micas = [[resnik(t1, t2)
                  for t2 in p2_terms]
                 for t1 in p1_terms]

        row_max = [max(row) for row in micas]
        col_max = [max([row[i] for row in micas]) for i in range(len(micas[0]))]

        if not scores or 'resnik' in scores:
            # average, max, best-match-average
            out['resnik_avg'] = sum([sum(row) for row in micas]) / (len(micas) * len(micas[0]))
            out['resnik_best_avg'] = sum(row_max) / len(row_max)

        if not scores or 'owlsim' in scores:
            owl_max_score = max(row_max)
            owl_avg_score = (sum(row_max) + sum(col_max)) / (len(row_max) + len(col_max))
            opt_p1_scores = [hpoic.get_term_ic(t) for t in p1_terms]
            opt_max_score = max(opt_p1_scores)
            opt_avg_score = sum(opt_p1_scores) / len(opt_p1_scores)

            out['owlsim_max'] = 100 * owl_max_score / opt_max_score
            out['owlsim_avg'] = 100 * owl_avg_score / opt_avg_score
            out['owlsim_combined'] = 0.5 * (out['owlsim_max'] + out['owlsim_avg'])

    if not scores or 'lin' in scores:
        lins = [[lin(t1, t2, res)
                 for (t2, res) in zip(p2_terms, row)]
                for (t1, row) in zip(p1_terms, micas)]

        out['lin_avg'] = sum([sum(row) for row in lins]) / (len(lins) * len(lins[0]))
        out['lin_best_avg'] = sum([max(row) for row in lins]) / len(lins)

    if not scores or 'jc' in scores:
        jcs = [[jc(t1, t2, res)
                for (t2, res) in zip(p2_terms, row)]
               for (t1, row) in zip(p1_terms, micas)]
        out['jc_avg'] = sum([sum(row) for row in jcs]) / (len(jcs) * len(jcs[0]))
        out['jc_best_avg'] = sum([max(row) for row in jcs]) / len(jcs)

    return out

def script(patient_filename, hpo_filename, disease_phenotype_filename,
           orphanet_lookup_filename=None, orphanet_prevalence_filename=None, proto=None,
           use_disease_prevalence=False, use_phenotype_frequency=False,
           use_patient_phenotypes=False, distribute_ic_to_leaves=False,
           use_external_ids=False, use_aoo=False, scores=None):
    hpo = HPO(hpo_filename, new_root='HP:0000118')
    diseases = Diseases(disease_phenotype_filename)

    orphanet = None
    if orphanet_lookup_filename and orphanet_prevalence_filename:
        orphanet = Orphanet(orphanet_prevalence_filename, lookup_filename=orphanet_lookup_filename)

    patients = [patient
                for patient in Patient.iter_from_file(patient_filename, hpo)
                if patient.hp_terms]

    if proto:
        proto = [patient
                 for patient in Patient.iter_from_file(proto, hpo)
                 if patient.hp_terms]

    if use_patient_phenotypes:
        use_patient_phenotypes = patients

    hpoic = HPOIC(hpo, diseases, orphanet=orphanet, patients=use_patient_phenotypes,
                  use_disease_prevalence=use_disease_prevalence,
                  use_phenotype_frequency=use_phenotype_frequency,
                  distribute_ic_to_leaves=distribute_ic_to_leaves)

    total_patient_logprob = 0
    for patient in patients:
        total_patient_logprob += hpoic.information_content(patient.hp_terms)

    logging.info('Total patient logprob: {:.1f}'.format(-total_patient_logprob))

    def get_id(patient):
        return patient.external_id if (patient.external_id and use_external_ids) else patient.id

    header = None
    for i in range(len(patients)):
        patient = patients[i]
        id1 = get_id(patient)
        compare_against = [patients[j] for j in range(i+1, len(patients))]
        if proto:
            compare_against.extend(proto)

        for o in compare_against:
            id2 = get_id(o)
            sims = compare_patients(hpoic, patient, o, scores=scores, use_aoo=use_aoo)
            if header is None:
                header = sorted(sims)
                print('\t'.join(['A', 'B'] + header))

            sim_strs = ['{:.6f}'.format(sims[sim]) for sim in header]
            for sim, sim_str in zip(header, sim_strs):
                logging.debug('{}: {}'.format(sim, sim_str))
            print('\t'.join([id1, id2] + sim_strs))


def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()

    parser = ArgumentParser(description=description)
    parser.add_argument('patient_filename', metavar='patients.json')
    parser.add_argument('hpo_filename', metavar='hp.obo')
    parser.add_argument('disease_phenotype_filename', metavar='phenotype_annotations.tab')
    parser.add_argument('--orphanet-lookup', metavar='en_product1.xml',
                        dest='orphanet_lookup_filename', default=None)
    parser.add_argument('--orphanet-prevalence', metavar='en_product2.xml',
                        dest='orphanet_prevalence_filename', default=None)
    parser.add_argument('--use-disease-prevalence', default=False, action='store_true')
    parser.add_argument('--use-phenotype-frequency', default=False, action='store_true')
    parser.add_argument('--use-patient-phenotypes', default=False, action='store_true')
    parser.add_argument('--distribute-ic-to-leaves', default=False, action='store_true')
    parser.add_argument('--use-aoo', default=False, action='store_true')
    parser.add_argument('--use-external-ids', default=False, action='store_true',
                        help='Use external ids, where available, instead of report ids')
    parser.add_argument('--log', dest='loglevel', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default='WARNING')
    parser.add_argument('--proto', metavar="file", help="HPO file of disease prototypes to compare against as well")
    parser.add_argument('-s', '--score', dest='scores', action='append', default=[],
                        choices=['jaccard', 'resnik', 'lin', 'jc', 'owlsim', 'jz', 'ui', 'simgic', 'nsimgic', 'icca'],
                        help='Include this score in the output for each pair of patients')

    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = parse_args(args)
    logging.basicConfig(level=args.__dict__.pop('loglevel'))
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
