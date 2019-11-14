#!/usr/bin/env python3

"""
Module to provide a class-based interface for interacting with diseases
and associated phenotypes (from a phenotype_annotations.tab file)
"""

__author__ = 'Orion Buske'

import os
import sys
import re
import logging

from collections import defaultdict

DBS = set(['ORPHANET', 'OMIM', 'DECIPHER'])
db_re = re.compile('([A-Z]+:\d+)')
FREQUENCIES = {'very rare':  0.01,
               'rare':       0.05,
               'occasional': 0.075,
               'frequent':   0.33,
               'typical':    0.5,
               'variable':   0.5,
               'common':     0.75,
               'hallmark':   0.9,
               'obligate':   1.0}
FREQUENCY_TERMS = {
    'HP:0040280': 1.0,  # Obligate
    'HP:0040281': (0.99 + 0.80) / 2.0,  # Very frequent
    'HP:0040282': (0.79 + 0.30) / 2.0,  # Frequent
    'HP:0040283': (0.05 + 0.29) / 2.0,  # Occasional
    'HP:0040284': (0.01 + 0.04) / 2.0,  # Very rare
    'HP:0040285': 0.0,  # Excluded
}
fraction_frequency_re = re.compile(r'of|/')

logger = logging.getLogger(__name__)


class Disease:
    def __init__(self, db, id, phenotype_freqs):
        self.db = db
        self.id = id
        self.phenotype_freqs = phenotype_freqs

    def __str__(self):
        return '{}:{}'.format(self.db, self.id)


class Diseases:
    def __init__(self, filename, db=None):
        """Read disease informaton from file into dict: {(db, id) -> Disease}

        db: database to restrict to (e.g. 'OMIM'), or None to include all

        Will use disease links in reference column (#6)
        """
        if db:
            assert db in DBS
            logger.info('Filtering to diseases in {}'.format(db))

        diseases = dict(self.iter_diseases(filename, filter_db=db))

        for i, disease in zip(range(5), diseases):
            logger.debug(diseases[disease].__dict__)

        self.diseases = diseases

    def __iter__(self):
        return iter(self.diseases)

    def __len__(self):
        return len(self.diseases)

    @classmethod
    def parse_frequency(cls, s, default=None):
        """Return float parsed frequency or default if problem or absent"""
        if s.upper() in FREQUENCY_TERMS:
            return FREQUENCY_TERMS[s]

        s = s.lower()
        if not s:
            freq = default
        elif s in FREQUENCIES:
            freq = FREQUENCIES[s]
        elif s.endswith('%'):
            s = s.replace('%', '')
            if '-' in s:
                # Average any frequency ranges
                low, high = s.split('-')
                freq = (float(low) + float(high)) / 2 / 100
            else:
                freq = float(s) / 100
        else:
            try:
                num, denom = fraction_frequency_re.split(s)
            except:
                logger.error("Error parsing frequency: {!r}".format(s))
                freq = default
            else:
                freq = float(num) / float(denom)

        return freq

    @classmethod
    def iter_diseases(cls, filename, default_freq=None, filter_db=None):
        disease_phenotypes = defaultdict(dict)  # (db, id) -> {HP: freq, ...}
        with open(filename, encoding='utf-8') as ifp:
            for line in ifp:
                line = line.rstrip()
                tokens = line.split('\t')
                if len(tokens) == 1: continue
                diseases = [(tokens[0].strip(), tokens[1].strip())]
                alt_diseases = db_re.findall(tokens[5].strip())
                #name = tokens[2].strip().split(';')[0]
                for alt_disease in alt_diseases:
                    db, id = alt_disease.split(':')
                    db = db.strip()
                    if db in set(['IM', 'MIM']):
                        db = 'OMIM'

                    id = int(id.strip())
                    if db not in DBS:
                        continue

                    diseases.append((db, id))

                hp_term = tokens[4].strip()
                freq = cls.parse_frequency(tokens[8])
                for disease in diseases:
                    if filter_db is not None and disease[0] != filter_db: continue

                    phenotypes = disease_phenotypes[disease]
                    if freq is not None and hp_term in phenotypes:
                        old_freq = phenotypes[hp_term]
                        # Always take max annotated frequency
                        if old_freq is None or old_freq < freq:
                            phenotypes[hp_term] = freq

                        if old_freq != freq:
                            logging.warning('Found conflicting frequencies ({}, {}) for same disease-phenotype: {}:{} - {} (taking larger)'.format(old_freq, freq, disease[0], disease[1], hp_term))
                    else:
                        phenotypes[hp_term] = freq

        for (db, id), phenotypes in disease_phenotypes.items():
            disease = Disease(db, id, phenotypes)
            yield ((db, id), disease)


def script(phenotype_filename):
    diseases = Diseases(phenotype_filename, db='ORPHANET')

    for id, disease in diseases.diseases.items():
        freqs = list(disease.phenotype_freqs.values())
        numeric = [f for f in freqs if f is not None]
        if len(numeric) >= 5 and max(numeric) < 0.1:
            logging.warning('Disease only has rare phenotypes with frequencies {}:{}'.format(id[0], id[1]))


def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()

    parser = ArgumentParser(description=description)
    parser.add_argument('phenotype_filename', metavar='phenotype_annotations.tab')

    return parser.parse_args(args)


def main(args=sys.argv[1:]):
    logging.basicConfig(level='INFO')
    args = parse_args(args)
    script(**vars(args))


if __name__ == '__main__':
    sys.exit(main())
