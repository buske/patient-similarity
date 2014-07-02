#!/usr/bin/env python

"""

"""

__author__ = 'Orion Buske'

import os
import sys
import logging

from math import log, exp
from collections import defaultdict
from collections import Counter

from hpo import HPO
from mim import MIM
from orphanet import Orphanet

EPS = 1e-9

def load_hpo(hpo_filename):
    hpo = HPO(hpo_filename)
    #hpo.filter_to_descendants('HP:0000001')
    hpo.filter_to_descendants('HP:0000118')
    logging.info('Found {} terms'.format(len(hpo)))
    return hpo

def _bound(p, eps=EPS):
    return min(max(p, eps), 1-eps)

class HPOIC:
    def __init__(self, hpo, mim, orphanet, patients=None,
                 use_disease_prevalence=False,
                 use_phenotype_frequency=False):
        logging.info('HPO root: {}'.format(hpo.root.id))

        term_freq = self.get_term_frequencies(mim, hpo, orphanet, patients=patients,
                                              use_disease_prevalence=use_disease_prevalence, 
                                              use_phenotype_frequency=use_phenotype_frequency)
        logging.info('Total term frequency mass: {}'.format(sum(term_freq.values())))

        term_ic = self.get_ics(hpo.root, term_freq)
        logging.info('IC calculated for {}/{} terms'.format(len(term_ic), len(hpo)))
        for term, ic in term_ic.items():
            for p in term.parents:
                assert ic >= term_ic[p] - EPS, str((term, ic, term_ic[p], p))

        lss = self.get_link_strengths(hpo.root, term_ic)
        logging.info('Link strength calculated for {}/{} terms'.format(len(lss), len(hpo)))

        t = hpo['HP:0000291']
        logging.info('{}: {}, {}'.format(t, term_ic[t], sum([lss[x] for x in t.ancestors() - set([t])])))
        for p in t.parents:
            logging.info('  {}: {}'.format(p, term_ic[p]))
            for c in p.children:
                logging.info('    {}: {}'.format(c, term_ic[c]))

        self.term_ic = term_ic
        self.lss = lss

#         for term in hpo:
#             ic = self.get_term_ic(term)
#             ls_ic = sum([lss.get(a, 0) for a in term.ancestors()])
#             logging.info('ic: {:.6f}, lsic: {:.6f}, term: {}'.format(ic, ls_ic, term))
#             assert abs(ic - ls_ic) < EPS * 2


    @classmethod
    def get_term_frequencies(cls, mim, hpo, orphanet=None, patients=None,
                             use_disease_prevalence=False, 
                             use_phenotype_frequency=False):
        bound = _bound
        raw_freq = defaultdict(float)

        if use_phenotype_frequency:
            # Use average observed phenotype frequency as default
            default_hp_freq = cls.get_average_phenotype_frequency(mim, hpo)
            logging.info('Average observed phenotype frequency: {:.4f}'.format(default_hp_freq))

        if use_disease_prevalence:
            default_disease_freq = orphanet.average_frequency()
            logging.info('Average disease frequency: {}'.format(default_disease_freq))

        # Aggregate weighted term frequencies across OMIM for IC corpus
        for disease in mim:
            # Weight by disease prevalence?
            if use_disease_prevalence:
                prevalence = orphanet.prevalence.get(disease.id)
                if prevalence is None:
                    prevalence = default_disease_freq
            else:
                prevalence = 1

            for hp_term, freq in disease.phenotype_freqs.items():
                # Try to resolve term
                try:
                    term = hpo[hp_term]
                except KeyError:
                    continue

                # Weight by phenotype frequency?
                if use_phenotype_frequency:
                    if freq is None:
                        freq = default_hp_freq
                else:
                    freq = 1

                raw_freq[term] += freq * prevalence

#         logging.warn('DISTRIBUTING WEIGHT TO LEAVES')
#         def distribute_weight_towards_leaves(node):
#             freq = raw_freq[node] + 0.01
#             if node.children:
#                 # Evenly divide frequency amongst children
#                 logging.info('Distributing {:.4f} to children of {}'.format(freq, node))
#                 raw_freq[node] = 0
#                 freq_part = freq / len(node.children)
#                 for child in node.children:
#                     raw_freq[child] += freq_part
#                     distribute_weight_towards_leaves(child)

#         distribute_weight_towards_leaves(hpo.root)

        if patients:
            for patient in patients:
                for term in patient.hp_terms:
                    raw_freq[term] += 1

        # Normalize all frequencies to sum to 1
        term_freq = {}
        total_freq = sum(raw_freq.values())
        for term in raw_freq:
            assert term not in term_freq
            term_freq[term] = bound(raw_freq[term] / total_freq)

        return term_freq

    @classmethod
    def get_descendant_lookup(cls, root, accum=None):
        if accum is None: accum = {}
        if root in accum: return

        descendants = set([root])
        for child in root.children:
            cls.get_descendant_lookup(child, accum)
            descendants.update(accum[child])
        accum[root] = descendants
        return accum

    @classmethod
    def get_ics(cls, root, term_freq):
        bound = _bound
        eps = EPS

        term_descendants = cls.get_descendant_lookup(root)

        term_ic = {}
        for node, descendants in term_descendants.items():
            prob_mass = 0.0
            for descendant in descendants:
                prob_mass += term_freq.get(descendant, 0)

            if prob_mass > eps:
                prob_mass = bound(prob_mass)
                term_ic[node] = -log(prob_mass)

        return term_ic
    
    @classmethod
    def get_link_strengths(cls, root, term_ic):
        # Calculate closer link strength (conditional IC) approximation to handle
        # multiple parents per node
        eps = EPS
        lss = {}
        for term in term_ic:
            assert term not in lss
            ls = 0
            if term.parents:
                ic = term_ic[term]
                ls = ic - max([term_ic[p] for p in term.parents])

            lss[term] = ls
            assert ls >= -eps, ls

#         def recurse(node):
#             if node in term_ic:
#                 ancestor_ic = 0
#                 for a in node.ancestors():
#                     if a != node:
#                         try:
#                             ancestor_ic += lss[a]
#                         except KeyError:
#                             return

#                 lss[node] = term_ic[node] - ancestor_ic
#                 for child in node.children:
#                     recurse(child)

#         recurse(root)

#         for term, ic in term_ic.items():
#             assert term not in lss
#             if len(term.parents) == 1:
#                 ls = term_ic[next(iter(term.parents))] - ic

            
#             siblings = None
#             for p in term.parents:
#                 if siblings is None:
#                     siblings = set(p.children)
#                 else:
#                     siblings.intersection_update(p.children)

#             if siblings is None:
#                 lss[term] = 0.0
#                 assert term.is_root(), 'Missing siblings for term: {!r}'.format(term)
#             else:
#                 # Siblings now contains all other siblings to term (+ self)
#                 sibling_prob_mass = 0.0
#                 for t in siblings:
#                     if t in term_ic:
#                         sibling_prob_mass += exp(-term_ic[t])

#                 if sibling_prob_mass < eps:
#                     lss[term] = 0.0
#                 else:
#                     lss[term] = term_ic[term] + log(sibling_prob_mass)
#                     assert lss[term] > 0

        return lss

    @classmethod
    def get_average_phenotype_frequency(cls, mim, hpo):
        freq_sum = 0
        n_freqs = 0
        dropped = set()
        for disease in mim:
            for hp_term, freq in disease.phenotype_freqs.items():
                try:
                    hpo[hp_term]
                except KeyError:
                    dropped.add(hp_term)
                else:
                    if freq:
                        freq_sum += freq
                        n_freqs += 1

        if dropped:
            logging.warning('Dropped {} unknown terms when computing average frequency'.format(len(dropped)))
        return freq_sum / n_freqs

    def get_term_ic(self, term):
        """Return information content of given term, falling back to parents as necessary"""
        ic = self.term_ic.get(term)
        if ic is None:
            if term.parents:
                ic = max([self.get_term_ic(p) for p in term.parents])
            else:
                ic = 0.0
        return ic

    def information_content(self, terms):
        """Return the information content of the given terms"""
        return sum([self.get_term_ic(term) for term in terms])

    def ls_information_content(self, ancestors):
        """Return the "joint" information content of the given bag of terms"""
        return sum([self.lss.get(term, 0) for term in ancestors])

    def counter_ls_information_content(self, ancestors):
        """Return the "joint" information content of the given counter of terms"""
        return sum([self.lss.get(term, 0) * count for term, count in ancestors.items()])


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

def compare_patients(hpoic, patient1, patient2, scores=None):
    logging.debug('Comparing patients: {}, {}'.format(patient1.id, patient2.id))

    p1_terms = patient1.hp_terms
    p2_terms = patient2.hp_terms
    assert p1_terms and p2_terms
    out = {}

#     logging.debug('Patient 1 terms and IC')
#     for t in p1_terms:
#         logging.debug('  {:.6f}: {}'.format(hpoic.get_term_ic(t), t))

#     logging.debug('Patient 2 terms and IC')
#     for t in p2_terms:
#         logging.debug('  {:.6f}: {}'.format(hpoic.get_term_ic(t), t))

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
        out['simgic_ls'] = hpoic.ls_information_content(common_ancestors) / hpoic.ls_information_content(all_ancestors)

    if not scores or 'jz' in scores:
        p1_ancestor_counts = ancestor_counts(p1_terms)
        p2_ancestor_counts = ancestor_counts(p2_terms)
        common_ancestor_counts = p1_ancestor_counts & p2_ancestor_counts  # min
        all_ancestor_counts = p1_ancestor_counts | p2_ancestor_counts  # max
        out['simgic_jz'] = 2 * hpoic.counter_ls_information_content(common_ancestor_counts) / (hpoic.information_content(p1_terms) + hpoic.information_content(p2_terms))
        out['simgic_jz2'] = hpoic.counter_ls_information_content(common_ancestor_counts) / hpoic.counter_ls_information_content(all_ancestor_counts)
        
    if not scores or 'jaccard' in scores:
        jaccard_rows = [[jaccard(t1, t2) for t2 in p2_terms] for t1 in p1_terms]
        out['jaccard_best_avg'] = sum([max(row) for row in jaccard_rows]) / len(jaccard_rows)
        out['jaccard_avg'] = sum([sum(row) for row in jaccard_rows]) / (len(jaccard_rows) * len(jaccard_rows[0]))
        out['jaccard_max'] = max([max(row) for row in jaccard_rows])


    if not scores or ('resnik' in scores or 'owlsim' in scores):
        micas = [[resnik(t1, t2) 
                  for t2 in p2_terms]
                 for t1 in p1_terms]

        row_max = [max(row) for row in micas]
        col_max = [max([row[i] for row in micas]) for i in range(len(micas[0]))]

        if 'resnik' in scores:
            # average, max, best-match-average
            out['resnik_avg'] = sum([sum(row) for row in micas]) / (len(micas) * len(micas[0]))
            out['resnik_best_avg'] = sum(row_max) / len(row_max)
            out['resnik_max'] = max(row_max)
            
        if 'owlsim' in scores:
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
        out['lin_max'] = max([max(row) for row in lins])

    if not scores or 'jc' in scores:
        jcs = [[jc(t1, t2, res)
                for (t2, res) in zip(p2_terms, row)]
               for (t1, row) in zip(p1_terms, micas)]
        out['jc_avg'] = sum([sum(row) for row in jcs]) / (len(jcs) * len(jcs[0]))
        out['jc_best_avg'] = sum([max(row) for row in jcs]) / len(jcs)
        out['jc_max'] = max([max(row) for row in jcs])

    if not scores or 'ob' in scores:
        p1_ic = hpoic.ls_information_content(p1_ancestors)
        p2_ic = hpoic.ls_information_content(p2_ancestors)
        shared_ic = hpoic.ls_information_content(common_ancestors)
        out['ob'] = 2 * shared_ic / (p1_ic + p2_ic)  # harmonic mean
        out['ob2'] = shared_ic / hpoic.ls_information_content(all_ancestors)

        logging.debug('Patient 1 ic: {:.6f}'.format(p1_ic))
        logging.debug('Patient 2 ic: {:.6f}'.format(p2_ic))
        logging.debug('Shared ic: {:.6f}'.format(shared_ic))

#     for t1 in p1_terms:
#         res, t2 = max([(resnik(t1, t2), t2) for t2 in p2_terms])
#         logging.debug('{} vs (best) {}'.format(t1, t2))
#         logging.debug('  mica: {}'.format(mica(t1, t2))) 
#         logging.debug('  mica: {}'.format(mica(t1, t2))) 
#         logging.debug('  resnik: {:.4f}'.format(resnik(t1, t2))) 
#         logging.debug('  lin: {:.4f}'.format(lin(t1, t2))) 
#         logging.debug('  j&c: {:.4f}'.format(jc(t1, t2))) 
#         logging.debug('  jaccard: {:.4f}'.format(jaccard(t1, t2))) 

    return out
        
class Patient:
    def __init__(self, id, hp_terms, neg_hp_terms=None, onset=None):
        self.id = id
        self.hp_terms = hp_terms
        self.neg_hp_terms = neg_hp_terms
        self.onset = onset
        self._ancestors = None

    def __repr__(self):
        return self.id

    def __lt__(self, o):
        return self.id < o.id

    def ancestors(self):
        if self._ancestors is None:
            ancestors = set()
            for term in self.hp_terms:
                ancestors.update(term.ancestors())
            self._ancestors = ancestors

        return self._ancestors

    @classmethod
    def from_line(cls, line, hpo, missing_terms=None):
        def resolve_terms(terms, missing_terms=missing_terms):
            nodes = []
            for term in terms:
                term = term.strip()
                try:
                    node = hpo[term]
                except KeyError:
                    if missing_terms is not None:
                        missing_terms.add(term)
                else:
                    nodes.append(node)
            return nodes

        entry = dict(zip(['id', 'hps', 'no_hps', 'onset'], line.strip().split('\t')))
        id = entry['id']
        hp_terms = entry.get('hps', [])
        if hp_terms:
            hp_terms = resolve_terms(hp_terms.split(','))

        neg_hp_terms = entry.get('no_hps', [])
        if neg_hp_terms:
            neg_hp_terms = resolve_terms(neg_hp_terms.split(','))

        onset = entry.get('onset')
        if onset:
            onset = resolve_terms(onset.split(','))

        return Patient(id, hp_terms, neg_hp_terms)
        

    @classmethod
    def iter_from_file(cls, filename, hpo):
        missing_terms = set()

        with open(filename) as ifp:
            for line in ifp:
                yield cls.from_line(line, hpo, missing_terms)

        if missing_terms:
            logging.warning('Could not find {} terms: {}'.format(len(missing_terms), ','.join(missing_terms)))


def script(patient_hpo_filename, hpo_filename, disease_phenotype_filename, 
           orphanet_lookup_filename, orphanet_prevalence_filename, proto=None, 
           use_disease_prevalence=False, use_phenotype_frequency=False, 
           use_patient_phenotypes=False, scores=None):
    hpo = load_hpo(hpo_filename)
    mim = MIM(disease_phenotype_filename)
    orphanet = Orphanet(orphanet_lookup_filename, orphanet_prevalence_filename)

    patients = [patient 
                for patient in Patient.iter_from_file(patient_hpo_filename, hpo)
                if patient.hp_terms]

    if proto:
        proto = [patient 
                 for patient in Patient.iter_from_file(proto, hpo)
                 if patient.hp_terms]

    if use_patient_phenotypes:
        use_patient_phenotypes = patients

    hpoic = HPOIC(hpo, mim, orphanet, patients=use_patient_phenotypes,
                  use_disease_prevalence=use_disease_prevalence,
                  use_phenotype_frequency=use_phenotype_frequency)

    total_patient_logprob = 0
    total_patient_ls_logprob = 0
    for patient in patients:
        total_patient_logprob += hpoic.information_content(patient.hp_terms)
        total_patient_ls_logprob += hpoic.ls_information_content(patient.ancestors())

    logging.info('Total patient logprob: {:.1f}'.format(-total_patient_logprob))
    logging.info('Total patient ls logprob: {:.1f}'.format(-total_patient_ls_logprob))

    header = None
    for i in range(len(patients)):
        patient = patients[i]
        compare_against = [patients[j] for j in range(i+1, len(patients))]
        if proto:
            compare_against.extend(proto)

        for o in compare_against:
            sims = compare_patients(hpoic, patient, o, scores=scores)
            if header is None:
                header = sorted(sims)
                print('#A\tB\t{}'.format('\t'.join(header)))

            sim_strs = ['{:.6f}'.format(sims[sim]) for sim in header]
            for sim, sim_str in zip(header, sim_strs):
                logging.debug('{}: {}'.format(sim, sim_str))
            print('\t'.join([patient.id, o.id] + sim_strs))


def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()
    
    parser = ArgumentParser(description=description)
    parser.add_argument('patient_hpo_filename', metavar='patients.hpo')
    parser.add_argument('hpo_filename', metavar='hp.obo')
    parser.add_argument('disease_phenotype_filename', metavar='phenotype_annotations.tab')
    parser.add_argument('orphanet_lookup_filename', metavar='orphanet_lookup')
    parser.add_argument('orphanet_prevalence_filename', metavar='orphanet_prevalence')
    parser.add_argument('--use-disease-prevalence', default=False, action='store_true')
    parser.add_argument('--use-phenotype-frequency', default=False, action='store_true')
    parser.add_argument('--use-patient-phenotypes', default=False, action='store_true')
    parser.add_argument('--log', dest='loglevel', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default='WARNING')
    parser.add_argument('--proto', metavar="file", help="Hpo file of disease prototypes to compare against as well")
    parser.add_argument('-s', '--score', dest='scores', action='append', default=[],
                        choices=['jaccard', 'resnik', 'lin', 'jc', 'owlsim', 'ob', 'jz', 'ui', 'simgic'],
                        help='Include this score in the output for each pair of patients')

    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = parse_args(args)
    logging.basicConfig(level=args.__dict__.pop('loglevel'))
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
