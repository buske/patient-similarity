#!/usr/bin/env python3

"""
"""


import sys
import os
import logging
import csv
import string

from collections import defaultdict

def read_cohorts(filename):
    cohort_lookup = {}
    with open(filename) as ifp:
        for line in ifp:
            line = line.strip()
            cohorts = set(line.strip().split())
            for cohort in cohorts:
                if cohort in cohort_lookup:
                    cohort_lookup[cohort].update(cohorts)
                else:
                    cohort_lookup[cohort] = cohorts
    return cohort_lookup

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

def is_same_cohort(c1, c2, cohort_lookup):
    return c1 == c2 or (c1 in cohort_lookup and c2 in cohort_lookup[c1])

def script(pheno_sim, cohort_file, id_file=None):
    cohort_lookup = read_cohorts(cohort_file)

    ids = None
    if id_file:
        ids = read_id_lookup(id_file)

    seen_ids = set()
    matches_found = defaultdict(set)
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
            same_cohort = is_same_cohort(c1, c2, cohort_lookup)

            seen_ids.add(p1)
            seen_ids.add(p2)
            if same_cohort:
                matches_found[p1].add(p2)
                matches_found[p2].add(p1)

            print('{:d}'.format(same_cohort))

    for id in sorted(seen_ids):
        if id in matches_found:
            logging.info('Found matches for {}: {}'.format(id, ','.join(sorted(matches_found[id]))))
        else:
            logging.info('No matches found for {}'.format(id))

    logging.info('Found {} IDs'.format(len(seen_ids)))
    logging.info('{} with matches'.format(len(matches_found)))
    logging.info('{} without matches'.format(len(seen_ids) - len(matches_found)))

def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()

    parser = ArgumentParser(description=description)
    parser.add_argument("--id-file")
    parser.add_argument("pheno_sim")
    parser.add_argument("cohort_file", metavar="cohorts.txt",
                        help="One cohort per line, whitespace separated patients")

    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()
