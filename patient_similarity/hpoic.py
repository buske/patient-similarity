#!/usr/bin/env python3

"""
Module that provides information-content functionality for the HPO.
"""

__author__ = 'Orion Buske (buske@cs.toronto.edu)'

import sys
import logging

from math import log, exp
from collections import defaultdict


EPS = 1e-9
logger = logging.getLogger(__name__)

def _bound(p, eps=EPS):
    return min(max(p, eps), 1-eps)

class HPOIC:
    def __init__(self, hpo, diseases, orphanet=None, patients=None,
                 use_disease_prevalence=False,
                 use_phenotype_frequency=False,
                 distribute_ic_to_leaves=False):
        logger.info('HPO root: {}'.format(hpo.root.id))

        term_freq = self.get_term_frequencies(diseases, hpo, orphanet=orphanet, patients=patients,
                                              use_disease_prevalence=use_disease_prevalence,
                                              use_phenotype_frequency=use_phenotype_frequency,
                                              distribute_ic_to_leaves=distribute_ic_to_leaves)
        logger.info('Total term frequency mass: {}'.format(sum(term_freq.values())))

        term_ic = self.get_ics(hpo.root, term_freq)
        logger.info('IC calculated for {}/{} terms'.format(len(term_ic), len(hpo)))
        for term, ic in term_ic.items():
            for p in term.parents:
                assert ic >= term_ic[p] - EPS, str((term, ic, term_ic[p], p))

        lss = self.get_link_strengths(hpo.root, term_ic)
        logger.info('Link strength calculated for {}/{} terms'.format(len(lss), len(hpo)))

        self.term_ic = term_ic
        self.lss = lss

    @classmethod
    def get_term_frequencies(cls, diseases, hpo, orphanet=None, patients=None,
                             use_disease_prevalence=False,
                             use_phenotype_frequency=False,
                             distribute_ic_to_leaves=False):
        bound = _bound
        raw_freq = defaultdict(float)

        if use_phenotype_frequency:
            # Use average observed phenotype frequency as default
            default_hp_freq = cls.get_average_phenotype_frequency(diseases, hpo)
            logger.info('Average observed phenotype frequency: {:.4f}'.format(default_hp_freq))

        if use_disease_prevalence:
            default_disease_freq = orphanet.average_frequency()
            logger.info('Average disease frequency: {}'.format(default_disease_freq))

        # Aggregate weighted term frequencies across database for IC corpus
        n_phenotype_freqs_used = 0
        n_disease_prevalences_used = 0
        n_entries = 0
        for id, disease in diseases.diseases.items():
            # Weight by disease prevalence?
            if use_disease_prevalence:
                prevalence = orphanet.prevalence.get(disease.id)
                if prevalence is None:
                    prevalence = default_disease_freq
                else:
                    n_disease_prevalences_used += 1
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
                        n_phenotype_freqs_used += 1
                else:
                    freq = 1

                n_entries += 1
                raw_freq[term] += freq * prevalence

        if use_phenotype_frequency:
            logging.warning('Used phenotype frequencies for {}/{} disease-phenotype associations'.format(n_phenotype_freqs_used, n_entries))

        if use_disease_prevalence:
            logging.warning('Used disease prevalences for {}/{} diseases'.format(n_disease_prevalences_used, len(diseases)))

        if patients:
            for patient in patients:
                for term in patient.hp_terms:
                    raw_freq[term] += 1

        if distribute_ic_to_leaves:
            logger.warning('DISTRIBUTING WEIGHT TO LEAVES')
            q = [hpo.root]

            while q:
                node = q.pop()
                if node.children:
                    q.extend(node.children)
                    # Evenly divide frequency amongst children
                    freq = raw_freq[node]
                    logger.info('Distributing {:.4f} to children of {}'.format(freq, node))
                    raw_freq[node] = 0
                    freq_part = freq / len(node.children)
                    for child in node.children:
                        raw_freq[child] += freq_part

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
        eps = EPS
        lss = {}

        # P(term&parents) = P(term|parents) P(parents)
        # P(term|parents) = P(term&parents) / P(parents)
        # P(parents) = probmass(descendants)
        for term, ic in term_ic.items():
            assert term not in lss
            if term.parents:
                max_parent_ic = max([term_ic[parent] for parent in term.parents])
                ls = max(ic - max_parent_ic, 0.0)
            else:
                ls = 0.0

            lss[term] = ls

        return lss

    @classmethod
    def get_average_phenotype_frequency(cls, diseases, hpo):
        freq_sum = 0
        n_freqs = 0
        dropped = set()
        for id, disease in diseases.diseases.items():
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
            logger.warning('Dropped {} unknown terms when computing average frequency: {}'.format(len(dropped), ', '.join(sorted(dropped))))
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
        """Return the information content of the given terms, without backoff"""
        return sum([self.term_ic.get(term, 0) for term in terms])

    def ls_information_content(self, ancestors):
        """Return the "joint" information content of the given bag of terms"""
        return sum([self.lss.get(term, 0) for term in ancestors])

    def counter_ls_information_content(self, ancestors):
        """Return the "joint" information content of the given counter of terms"""
        return sum([self.lss.get(term, 0) * count for term, count in ancestors.items()])

    def counter_information_content(self, ancestors):
        """Return the sum of information content of the given counter of terms"""
        return sum([self.term_ic.get(term, 0) * count for term, count in ancestors.items()])
