#!/usr/bin/env python

"""
Representation of a patient, with phenotypes described with HPO terms.
"""

__author__ = 'Orion Buske (buske@cs.toronto.edu)'

import logging

logger = logging.getLogger(__name__)

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
            logger.warning('Could not find {} terms: {}'.format(len(missing_terms), ','.join(missing_terms)))
