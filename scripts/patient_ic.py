#!/usr/bin/env python3

"""
Helper script that prints out the information content for each patient
"""

__author__ = 'Orion Buske (buske@cs.toronto.edu)'

import sys
import logging

from hpo import HPO
from disease import Diseases
from hpoic import HPOIC
from patient import Patient

def script(patient_hpo_filename, hpo_filename, disease_phenotype_filename,
           **kwargs):
    hpo = HPO(hpo_filename, new_root='HP:0000118')
    diseases = Diseases(disease_phenotype_filename)
    hpoic = HPOIC(hpo, diseases)

    print('\t'.join(['Patient ID', 'External ID', 'IC']))
    for patient in Patient.iter_from_file(patient_hpo_filename, hpo):
        print('\t'.join(map(str, [patient.id, patient.external_id, hpoic.information_content(patient.hp_terms)])))

def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()

    parser = ArgumentParser(description=description)
    parser.add_argument('patient_hpo_filename', metavar='patients.hpo')
    parser.add_argument('hpo_filename', metavar='hp.obo')
    parser.add_argument('disease_phenotype_filename', metavar='phenotype_annotations.tab')
    parser.add_argument('--log', dest='loglevel', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default='WARNING')

    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = parse_args(args)
    logging.basicConfig(level=args.loglevel)

    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
