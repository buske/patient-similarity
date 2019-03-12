#!/usr/bin/env python3

"""
Representation of a patient, with phenotypes described with HPO terms.
"""

__author__ = 'Orion Buske (buske@cs.toronto.edu)'

import json
import logging
import csv

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
    def from_row(cls, row, hpo):
        pid = row['report_id']
        external_id = row.get('external_id')

        hpids = set()
        features = row.get('features', [])
        prenatal_perinatal_phenotype = row.get('prenatal_perinatal_phenotype', {})
        prenatal_phenotype = prenatal_perinatal_phenotype.get('prenatal_phenotype', [])
        features.extend(prenatal_phenotype)

        for feature in features:
            if feature.get('observed', 'yes') == 'yes':
                term = feature.get('id')
                if term:
                    hpids.add(term.strip())
                else:
                    logger.error('ID missing from term {}'.format(feature))

        for nonstandard in row.get('nonstandard_features', []):
            if nonstandard['observed'] == 'yes':
                for category in nonstandard.get('categories', []):
                    hpids.add(category['id'].strip())

        hp_terms = set()
        for hpid in hpids:
            try:
                term = hpo[hpid]
                hp_terms.add(term)
            except KeyError:
                logger.warning('Dropping term not found in HPO: {}'.format(hpid))

        # XXX: parse negative phenotypes

        onsets = []
        for term in row.get('global_age_of_onset', []):
            term_id = term['id']
            assert term_id.startswith('HP:') and len(term_id) == 10, 'Invalid onset: {!r}'.format(term_id)
            onsets.append(term_id)

        if len(onsets) > 1:
            logger.warning('Found multiple ages of onset... choosing randomly')

        onset = None
        if onsets:
            onset = onsets[0]

        diagnoses = set()
        for diagnosis in row.get('disorders', []):
            term = diagnosis.get('id')
            if term:
                diagnoses.add(term.strip())

        return Patient(id=pid, hp_terms=hp_terms,
                       onset=onset, diagnoses=diagnoses, external_id=external_id)

    @classmethod
    def iter_from_file(cls, filename, hpo):
        with open(filename) as ifp:
            database = json.load(ifp)
            for row in database:
                yield cls.from_row(row, hpo)

    @classmethod
    def iter_from_csv_file(cls, filename, hpo):
        with open(filename, newline='') as ifp:
            reader = csv.reader(ifp)
            for row in reader:
                if len(row) < 2:
                    raise FormatError('Expected at least two columns: id, phenotype1, ...')
                pid, term_ids = [row[0], row[1:]]
                terms = set()
                for term_id in term_ids:
                    term_id = term_id.strip()
                    if not term_id:
                        continue

                    term = hpo[term_id]
                    if term:
                        terms.add(term)

                yield Patient(id=pid, hp_terms=terms)
