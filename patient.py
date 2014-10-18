#!/usr/bin/env python

"""
Representation of a patient, with phenotypes described with HPO terms.
"""

__author__ = 'Orion Buske (buske@cs.toronto.edu)'

import csv
import logging

logger = logging.getLogger(__name__)

class Patient:
    def __init__(self, id, hp_terms, neg_hp_terms=None, onset=None, diagnoses=None, external_id=None):
        self.id = id
        self.hp_terms = hp_terms
        self.neg_hp_terms = neg_hp_terms
        self.onset = onset
        self.diagnoses = diagnoses
        self.external_id = external_id
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
    def from_row(cls, row, hpo, missing_terms=None):
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

        pid = row['Patient ID']
        external_id = row['External ID']

        hp_terms = row['HPO+']
        if hp_terms:
            hp_terms = resolve_terms(hp_terms.split(','))

        neg_hp_terms = row['HPO-']
        if neg_hp_terms:
            neg_hp_terms = resolve_terms(neg_hp_terms.split(','))


        onset = row.get('AOO')
        if onset:
            assert onset.startswith('HP:') and len(onset) == 10, 'Invalid onset: {!r}'.format(onset)

        diagnoses = row.get('Diagnoses')
        if diagnoses:
            diagnoses = set(diagnoses.split(','))

        return Patient(id=pid, hp_terms=hp_terms, neg_hp_terms=neg_hp_terms, 
                       onset=onset, diagnoses=diagnoses, external_id=external_id)
        

    @classmethod
    def iter_from_file(cls, filename, hpo):
        missing_terms = set()

        with open(filename) as ifp:
            reader = csv.DictReader(ifp, delimiter='\t')
            for row in reader:
                yield cls.from_row(row, hpo, missing_terms)

        if missing_terms:
            logger.warning('Could not find {} terms: {}'.format(len(missing_terms), ','.join(missing_terms)))
