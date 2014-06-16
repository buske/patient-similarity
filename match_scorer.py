#!/usr/bin/env python

"""
Given all the information for a pair of patients, computer a match score, with gene and variant info
"""


import sys
import os
import math
import logging
import csv
from collections import defaultdict
from itertools import combinations, product
from numpy import array


solutions = """
173 493
248 86
200 502
490 498
408 468
178 233 C1038
258 444
440 C1039
171 330
C1011 C1020
109 242 381 C1026
242 84
UDP_2543 UDP_2542
"""

cohort_lookup = {}
for line in solutions.strip().split('\n'):
    cohorts = set(line.strip().split())
    for cohort in cohorts:
        if cohort in cohort_lookup:
            cohort_lookup[cohort].update(cohorts)
        else:
            cohort_lookup[cohort] = cohorts

def read_pc_ids(filename):
    ids = {}  # id -> external
    with open(filename) as ifp:
        reader = csv.DictReader(ifp)
        for row in reader:
            patient_id = row['Report ID'].strip()
            external_id = row['Identifier'].strip()
            if external_id:
                ids[patient_id] = external_id
            
    return ids

def get_cohort(p):
    return p if p.startswith('UDP') else p.split('_')[0]

def is_same_cohort(c1, c2):
    return c1 == c2 or (c1 in cohort_lookup and c2 in cohort_lookup[c1])

def script(pc_file, pheno_sim):
    ids = read_pc_ids(pc_file)

    with open(pheno_sim) as ifp:
        for line in ifp:
            if not line or line.startswith('#'): continue
            tokens = line.strip().split('\t')
            p1, p2 = tokens[:2]
            scores = list(map(float, tokens[2:]))

            p1 = ids.get(p1, p1)
            p2 = ids.get(p2, p2)
            c1 = get_cohort(p1)
            c2 = get_cohort(p2)
            print('{:d}\t{}\t{}\t{}'.format(is_same_cohort(c1, c2), p1, p2, '\t'.join(map(str, scores))))

def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()
    
    parser = ArgumentParser(description=description)
    parser.add_argument("pc_file")
    parser.add_argument("pheno_sim")

    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()
        


    
