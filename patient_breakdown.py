#!/usr/bin/env python3

"""

"""

__author__ = 'Orion Buske'

import os
import sys
import logging

from math import log, exp
from collections import defaultdict

from disease import Diseases
from orphanet import Orphanet
from patient_similarity import PatientComparator, Patient


def script(patient1, patient2, patient_hpo_filename, hpo_filename, 
           disease_phenotype_filename, 
           orphanet_lookup_filename, orphanet_prevalence_filename, 
           **kwargs):
    hpo = HPO(hpo_filename, new_root='HP:0000118')
    patients = dict([(patient.id, patient)
                     for patient in Patient.iter_from_file(patient_hpo_filename, hpo)
                     if patient.hp_terms])

    p1 = patients[patient1]
    p2 = patients[patient2]

    diseases = Diseases(disease_phenotype_filename)
    orphanet = Orphanet(orphanet_prevalence_filename, lookup_filename=orphanet_lookup_filename)
    comparator = PatientComparator(hpo, diseases, orphanet)
    
    clusters = comparator.similarity_breakdown(p1, p2)
    cluster_strs = []
    for cluster in clusters:
        score, root, p1_terms, p2_terms = cluster
        if p1_terms:
            p1_terms = '\n'.join(['\t{}: {}'.format(t.id, t.name) for t in p1_terms])
        else:
            p1_terms = '\tNo matching terms'

        if p2_terms:
            p2_terms = '\n'.join(['\t{}: {}'.format(t.id, t.name) for t in p2_terms])
        else:
            p2_terms = '\tNo matching terms'

        if root:
            root = '{}: {} ({:.2f})'.format(root.id, root.name, score)
        else:
            root = 'Unmatched terms'

        cluster_strs.append('{}\n\t-- {} --\n{}\n\t-- {} --\n{}'.format(root, patient1, p1_terms, patient2, p2_terms))

    print('\n\n'.join(cluster_strs))

def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()
    
    parser = ArgumentParser(description=description)
    parser.add_argument('patient1')
    parser.add_argument('patient2')
    parser.add_argument('patient_hpo_filename', metavar='patients.hpo')
    parser.add_argument('hpo_filename', metavar='hp.obo')
    parser.add_argument('disease_phenotype_filename', metavar='phenotype_annotations.tab')
    parser.add_argument('orphanet_lookup_filename', metavar='orphanet_lookup')
    parser.add_argument('orphanet_prevalence_filename', metavar='orphanet_prevalence')
    parser.add_argument('--log', dest='loglevel', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default='WARNING')

    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = parse_args(args)
    logging.basicConfig(level=args.loglevel)

    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
