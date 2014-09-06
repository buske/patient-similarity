#!/usr/bin/env python

"""
"""


import sys
import os
import math
import logging
import csv
import string

from collections import defaultdict
from itertools import combinations, product
from numpy import array


solutions = """
174 F0000009 F0000010 F0000011 F0000012 F0000013 F0000014 F0000015 F0000016 F0000017 F0000018 F0000019 F0000020 F0000021
MFDM C1038 233
UDP_1019 UDP_1995 UDP_3510 UDP_2058 UDP_4626
UDP_1248 UDP_887
UDP_1258 UDP_2467
UDP_1262 UDP_5630
UDP_1363 UDP_5808
UDP_1480 UDP_1481
UDP_2146 UDP_2156
UDP_2249 UDP_2248
UDP_2542 UDP_10008
UDP_2591 UDP_2590
UDP_2610 UDP_3306
UDP_2652 UDP_2653
UDP_2753 UDP_2755
UDP_2803 UDP_2805
UDP_2994 UDP_2996 UDP_1706
UDP_337 UDP_338
UDP_4106 UDP_2780
UDP_4382 UDP_3866
UDP_4399 UDP_5824
UDP_4574 UDP_4478
UDP_4753 UDP_3491
UDP_4771 UDP_4770
UDP_4817 UDP_3155
UDP_4945 UDP_357
UDP_5020 UDP_8279 UDP_5019
UDP_503 UDP_4927
UDP_5356 UDP_5357
UDP_5433 UDP_5434
UDP_5433 UDP_5736 UDP_5434
UDP_5454 UDP_4942
UDP_5493 UDP_5492 UDP_5075
UDP_5510 UDP_5228 UDP_5702
UDP_579 UDP_577
UDP_608 UDP_606
UDP_6120 UDP_7528
UDP_629 UDP_499 UDP_10086
UDP_6597 UDP_6599
UDP_6613 UDP_1008
UDP_675 UDP_2298
UDP_7097 UDP_7075
UDP_7107 UDP_7106
UDP_7230 UDP_7234
UDP_7552 UDP_7553 UDP_7556
UDP_761 UDP_760
UDP_7672 UDP_5940
UDP_8027 UDP_7537
UDP_809 UDP_810
UDP_8117 UDP_8116
UDP_930 UDP_929
"""

cohort_lookup = {}
for line in solutions.strip().split('\n'):
    cohorts = set(line.strip().split())
    for cohort in cohorts:
        if cohort in cohort_lookup:
            cohort_lookup[cohort].update(cohorts)
        else:
            cohort_lookup[cohort] = cohorts

def read_id_lookup(filename):
    ids = {}  # id -> external
    with open(filename) as ifp:
        reader = csv.DictReader(ifp, delimiter='\t')
        for row in reader:
            patient_id = row['Patient ID'].strip()
            external_id = row['External ID'].strip()
            if external_id:
                ids[patient_id] = external_id
            
    return ids

def get_cohort(p):
    if p.startswith('UDP') or p.startswith('DDN'):
        return p
    else:
        cohort = p.split('_')[0].rstrip(string.ascii_lowercase)
        if cohort == '455CA': return '455'
        return cohort

def is_same_cohort(c1, c2):
    return c1 == c2 or (c1 in cohort_lookup and c2 in cohort_lookup[c1])

def script(pheno_sim, id_file=None):
    ids = None
    if id_file:
        ids = read_id_lookup(id_file)

    with open(pheno_sim) as ifp:
        # Burn and replace header
        header = ifp.readline()
        print('resp')

        for line in ifp:
            line = line.strip()
            if not line: continue

            tokens = line.split('\t')
            p1, p2 = tokens[:2]
            
            if ids:
                p1 = ids.get(p1, p1)
                p2 = ids.get(p2, p2)

            c1 = get_cohort(p1)
            c2 = get_cohort(p2)
            
            print('{:d}'.format(is_same_cohort(c1, c2)))

def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()
    
    parser = ArgumentParser(description=description)
    parser.add_argument("--id-file")
    parser.add_argument("pheno_sim")

    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()
        


    
