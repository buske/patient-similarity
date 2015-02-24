#!/usr/bin/env python3

"""
"""


import sys
import os
import logging
import string
import js

from collections import defaultdict

def read_cohorts(filename):
    cohort_lookup = {}
    with open(filename) as ifp:
        for line in ifp:
            line = line.strip()
            # Two-step cohort merging
            cohort = set(line.strip().split())

            # Combine cohorts of all patients in cohort
            for patient in cohort:
                cohort.update(cohort_lookup.get(patient, set()))

            # Overwrite cohort for any patients
            for patient in cohort:
                cohort_lookup[patient] = cohort

    return cohort_lookup

def is_same_cohort(p1, p2, lookup):
    return p1 == p2 or (p1 in lookup and p2 in lookup and lookup[p1] == lookup[p2])

def script(matching_file, cohort_file):
    cohort_lookup = read_cohorts(cohort_file)

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

            same_cohort = is_same_cohort(p1, p2, cohort_lookup)

            seen_ids.add(p1)
            seen_ids.add(p2)
            if same_cohort:
                matches_found[p1].add(p2)
                matches_found[p2].add(p1)

            # Print boolean as binary
            print(int(same_cohort))

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
    parser.add_argument("matching_file", metavar="patients.sim")
    parser.add_argument("cohort_file", metavar="cohorts.txt",
                        help="One cohort per line, whitespace separated patient IDs")

    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()
