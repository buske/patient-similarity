#!/usr/bin/env python

"""

"""

__author__ = 'Orion Buske (buske@cs.toronto.edu)'

import sys
import logging
import pickle

from numpy import array
from random import sample

from hpo import HPO
from disease import Diseases
from patient import Patient
from hpoic import HPOIC
from patient_similarity import compare_patients

logger = logging.getLogger(__name__)
DEFAULT_TERMS = 1
DEFAULT_REPLICATES = 100000
DEFAULT_SCORE = 'icca'
TOP_FRAC = 0.1

def get_filename(filebase, n_terms):
    return '{}.{}.pkl'.format(filebase, n_terms)

def read_distribution(filename):
    with open(filename, 'rb') as ifp:
        return pickle.load(ifp)

def write_distribution(distribution, filename):
    with open(filename, 'wb') as ofp:
        pickle.dump(distribution, ofp)

def diseases_to_patients(diseases, hpo):
    disease_patients = []
    for disease in diseases:
        # Convert disease to patient
        disease_terms = []
        for term in disease.phenotype_freqs:
            try:
                term = hpo[term]
            except KeyError:
                continue

            disease_terms.append(term)

        if not disease_terms:
            logger.warn('Found no phenotypes for disease: {}'.format(disease))
            continue

        disease_patients.append(Patient(str(disease), disease_terms))

    return disease_patients

def calc_distribution(hpo, diseases, hpoic, n_replicates, n_terms, score):
    hpo_terms = set(hpo.hps.values())
    def get_random_patient(n):
        return Patient('random', sample(hpo_terms, n))

    # For each disease, an array distribution of similarity scores
    # Only the top 10% of n_replicates scores are stored
    distribution = {}  # disease -> array

    top_n = int(TOP_FRAC * n_replicates)
    logging.info('Creating {} patients with {} terms each...'.format(n_replicates, n_terms))
    patients = [get_random_patient(n_terms) for i in range(n_replicates)] 

    # Convert disease objects to patient objects for comparison
    diseases = diseases_to_patients(diseases, hpo)
    logging.info('Comparing against {} diseases...'.format(len(diseases)))

    for disease in diseases:
        sims = []
        for patient in patients:
            sim = compare_patients(hpoic, patient, disease, scores=[score])[score]
            sims.append(float(sim))

        sims.sort(reverse=True)
        top_scores = array(sims[:top_n])
        distribution[str(disease)] = top_scores

    return distribution

def script(hpo_filename, disease_phenotype_filename, 
           out_filebase, n_replicates=DEFAULT_REPLICATES, n_terms=DEFAULT_TERMS,
           score=DEFAULT_SCORE):
    hpo = HPO(hpo_filename, new_root='HP:0000118')
    diseases = Diseases(disease_phenotype_filename, db='ORPHANET')
    hpoic = HPOIC(hpo, diseases)

    distribution = calc_distribution(hpo, diseases, hpoic, n_replicates, n_terms, score)
    out_filename = get_filename(out_filebase, n_terms)
    write_distribution(distribution, out_filename)

def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()
    
    parser = ArgumentParser(description=description)
    parser.add_argument('hpo_filename', metavar='hp.obo')
    parser.add_argument('disease_phenotype_filename', metavar='phenotype_annotations.tab')
    parser.add_argument('out_filebase', metavar='out_filebase')
    parser.add_argument('--log', dest='loglevel', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default='WARNING')
    parser.add_argument('-t', dest='n_terms', default=DEFAULT_TERMS, type=int,
                        help='Number of terms to generate, per disease, per replicate')
    parser.add_argument('-n', dest='n_replicates', default=DEFAULT_REPLICATES, type=int,
                        help='Number of replicates to generate, per'
                        ' disease, per number of terms')
    parser.add_argument('-s', '--score', dest='score', default=DEFAULT_SCORE,
                        choices=['jaccard', 'resnik', 'lin', 'jc', 'owlsim', 'ob', 'jz', 'ui', 'simgic', 'icca'],
                        help='Generate a distribution for the specified score')

    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = parse_args(args)
    logging.basicConfig(level=args.__dict__.pop('loglevel'))
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
