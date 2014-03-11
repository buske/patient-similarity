#!/usr/bin/env python

"""

"""

__author__ = 'Orion Buske'

import os
import sys
import logging

from math import log, exp
from collections import defaultdict

from hpo import HPO
from mim import MIM
from orphanet import Orphanet

EPS = 1e-9


def load_hpo(hpo_filename):
    hpo = HPO(hpo_filename)
    hpo.filter_to_descendants('HP:0000118')
    #hpo.filter_to_descendants('HP:0000001')
    logging.info('Found {} terms'.format(len(hpo)))
    return hpo

class PatientComparator:
    def __init__(self, hpo, mim, orphanet):
        def bound(p, eps=EPS):
            return min(max(p, eps), 1-eps)

        raw_freq = defaultdict(float)
        freq_denom = 0

        # Use average observed phenotype frequency as default
        default_hp_freq = self.get_average_phenotype_frequency(mim, hpo)
        default_disease_freq = orphanet.average_frequency()
        logging.info('Average observed phenotype frequency: {:.4f}'.format(default_hp_freq))
        logging.info('Average disease frequency: {}'.format(default_disease_freq))
        for disease in mim:
            prevalence = orphanet.prevalence.get(disease.id)
            if prevalence is None:
                prevalence = default_disease_freq

            for hp_term, freq in disease.phenotype_freqs.items():
                try:
                    term = hpo[hp_term]
                except KeyError:
                    continue

                if freq is None:
                    freq = default_hp_freq

                weighted_freq = freq * prevalence
                freq_denom += weighted_freq
                raw_freq[term] += weighted_freq

        term_freq = {}
        for term in raw_freq:
            assert term not in term_freq
            term_freq[term] = bound(raw_freq[term] / freq_denom)

        def get_all_descendants(root, accum=None):
            if accum is None: accum = {}
            if root in accum: return

            descendants = set([root])
            for child in root.children:
                get_all_descendants(child, accum)
                descendants.update(accum[child])
            accum[root] = descendants
            return accum

        term_descendants = get_all_descendants(hpo.root)

        def get_ics(bound=bound, term_freq=term_freq, 
                    term_descendants=term_descendants):
            term_ic = {}
            for node, descendants in term_descendants.items():
                prob_mass = 0.0
                for descendant in descendants:
                    prob_mass += term_freq.get(descendant, 0)

                if prob_mass > EPS:
                    prob_mass = bound(prob_mass)
                    term_ic[node] = -log(prob_mass)

            return term_ic

        logging.info('HPO root: {}'.format(hpo.root.id))
        term_ic = get_ics()
        logging.info('IC calculated for {} terms'.format(len(term_ic)))

        # Calculate closer conditional IC approximation to handle multiple
        # parents per node
        ic_cond_parents = {}
        for term in term_ic:
            assert term not in ic_cond_parents
            siblings = None
            for p in term.parents:
                if siblings is None:
                    siblings = set(p.children)
                else:
                    siblings.intersection_update(p.children)

            if siblings is None:
                ic_cond_parents[term] = 0.0
                logging.info('Missing siblings for term: {!r}'.format(term))
            else:
                # Siblings now contains all other siblings to term
                sibling_prob_mass = 0.0
                for t in siblings:
                    if t in term_ic:
                        sibling_prob_mass += exp(-term_ic[t])

                if sibling_prob_mass < EPS:
                    ic_cond_parents[term] = 0.0
                else:
                    ic_cond_parents[term] = term_ic[term] + log(sibling_prob_mass)

        self.hpo = hpo
        self.mim = mim
        self.orphanet = orphanet
        self.term_freq = term_freq
        self.term_ic = term_ic
        self.ic_cond_parents = ic_cond_parents

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

    def bag_information_content(self, terms):
        """Return the information content of the given terms"""
        return sum([self.get_term_ic(term) for term in terms])

    def joint_information_content(self, ancestors):
        """Return the "joint" information content of the given bag of terms"""
        return sum([self.ic_cond_parents.get(term, 0) for term in ancestors])

    def compare(self, patient1, patient2):
        logging.debug('Comparing patients: {}, {}'.format(patient1.id, patient2.id))

        assert patient1.hp_terms and patient2.hp_terms

        logging.debug('Patient 1 terms and IC')
        for t in patient1.hp_terms:
            logging.debug('  {:.6f}: {} ({})'.format(self.get_term_ic(t), t, t.name))

        logging.debug('Patient 2 terms and IC')
        for t in patient2.hp_terms:
            logging.debug('  {:.6f}: {} ({})'.format(self.get_term_ic(t), t, t.name))

        p1_ancestors = patient1.ancestors()
        p2_ancestors = patient2.ancestors()
        common_ancestors = p1_ancestors & p2_ancestors  # min
        p1_ic = self.joint_information_content(p1_ancestors)
        p2_ic = self.joint_information_content(p2_ancestors)
        shared_ic = self.joint_information_content(common_ancestors)

        logging.debug('Patient 1 ic: {:.6f}'.format(p1_ic))
        logging.debug('Patient 2 ic: {:.6f}'.format(p2_ic))
        logging.debug('Shared ic: {:.6f}'.format(shared_ic))
        if shared_ic < EPS:
            harmonic_mean = 0.0
        else:
            harmonic_mean = 2 / (p1_ic / shared_ic + p2_ic / shared_ic)

        return [harmonic_mean, shared_ic]
        
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
    def iter_from_file(self, filename, hpo):
        missing_terms = set()
        def resolve_terms(terms, missing_terms=missing_terms):
            nodes = []
            for term in terms:
                term = term.strip()
                try:
                    node = hpo[term]
                except KeyError:
                    missing_terms.add(term)
                else:
                    nodes.append(node)
            return nodes

        with open(filename) as ifp:
            for line in ifp:
                entry = dict(zip(['id', 'hps', 'no_hps', 'onset'], line.strip().split('\t')))
                id = entry['id']
                hp_terms = entry.get('hps', [])
                if hp_terms:
                    hp_terms = resolve_terms(hp_terms.split(';'))

                neg_hp_terms = entry.get('no_hps', [])
                if neg_hp_terms:
                    neg_hp_terms = resolve_terms(neg_hp_terms.split(';'))
                    
                onset = entry.get('onset')
                if onset:
                    onset = resolve_terms(onset.split(';'))
                    
                yield Patient(id, hp_terms, neg_hp_terms)

        if missing_terms:
            logging.warning('Could not find {} terms: {}'.format(len(missing_terms), ','.join(missing_terms)))




def script(patient_hpo_filename, hpo_filename, disease_phenotype_filename, 
           orphanet_lookup_filename, orphanet_prevalence_filename, proto=None, 
           **kwargs):
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

    comparator = PatientComparator(hpo, mim, orphanet)

    scores = {}
    for i in range(len(patients)):
        for j in range(i+1, len(patients)):
            score = comparator.compare(patients[i], patients[j])
            scores[(i, j)] = score

    for i in range(len(patients)):
        for j in range(len(patients)):
            if i != j:
                score = scores[(min(i, j), max(i, j))]
                print('\t'.join([patients[i].id, patients[j].id, '{:.6f}\t{:.6f}'.format(score[0], score[1])]))

        if proto:
            for j in range(len(proto)):
                score = comparator.compare(patients[i], proto[j])
                print('\t'.join([patients[i].id, proto[j].id, '{:.6f}\t{:.6f}'.format(score[0], score[1])]))


def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()
    
    parser = ArgumentParser(description=description)
    parser.add_argument('patient_hpo_filename', metavar='patients.hpo')
    parser.add_argument('hpo_filename', metavar='hp.obo')
    parser.add_argument('disease_phenotype_filename', metavar='phenotype_annotations.tab')
    parser.add_argument('orphanet_lookup_filename', metavar='orphanet_lookup')
    parser.add_argument('orphanet_prevalence_filename', metavar='orphanet_prevalence')
    parser.add_argument('--log', dest='loglevel', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default='WARNING')
    parser.add_argument('--proto', metavar="file", help="hpo file of disease prototypes")

    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = parse_args(args)
    logging.basicConfig(level=args.loglevel)

    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
