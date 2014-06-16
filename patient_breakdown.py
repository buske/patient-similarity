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
from patient_similarity import PatientComparator, Patient, load_hpo


def script(patient_hpo_filename, hpo_filename, disease_phenotype_filename, 
           orphanet_lookup_filename, orphanet_prevalence_filename, proto=None, 
           **kwargs):
    hpo = load_hpo(hpo_filename)
    mim = MIM(disease_phenotype_filename)
    orphanet = Orphanet(orphanet_lookup_filename, orphanet_prevalence_filename)

    patients = dict([(patient.id, patient)
                     for patient in Patient.iter_from_file(patient_hpo_filename, hpo)
                     if patient.hp_terms])

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
                score_str = ['{:.6f}'.format(s) for s in score]
                print('\t'.join([patients[i].id, patients[j].id] + score_str))

        if proto:
            for j in range(len(proto)):
                score = comparator.compare(patients[i], proto[j])
                score_str = ['{:.6f}'.format(s) for s in score]
                print('\t'.join([patients[i].id, proto[j].id] + score_str))


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
