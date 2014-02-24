#!/usr/bin/env python

import sys
import numpy as np
import argparse
import pprint

parser = argparse.ArgumentParser(description='Reads in a list of patients. Each line of input should start with a patient identifier and then a series of HPO terms separated by semicolons describing that patient. The algorithm computes the similarity between all pairs of patients. It takes a bit of time to load up some pre-computed values.')
parser.add_argument('-x', dest='x', action='store_true', help='If this flag is on, the algorithm uses a slightly better bound on the conditional probabilities. See the readme for details.')
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='Input file. If not provided, reads from stdin.')
args = parser.parse_args()

from all_data import *
from phenotype import *
import sim


def load_patient(s):
    codes = []
    for hpid in s.strip().split(';'):
        if hpid in code_to_term:
            codes.append(code_to_term[hpid])
        else:
            print("Skipping unknown term: %s" % hpid, file=sys.stderr)

    return Patient(codes)


# Read in data
ifp = args.infile
patients = []
pids = []
for line in ifp:
    line = line.strip()
    if line.startswith('#'): continue

    if line:
        pid, codes = line.split(None, 1)
        pids.append(pid)
        patients.append(load_patient(codes))


sims = {}
# Compute pairwise scores
for i in range(len(patients) - 1):
    for j in range(i + 1, len(patients)):
        sims[(i, j)] = sim.sim(patients[i], patients[j], args.x)[0]
        # Print out pairwise results
        score = sims[(i, j)]
        print('\t'.join([pids[i], pids[j], '{:.6f}'.format(score)]))

# # Print out results
# print('\t'.join(pids))
# for i in range(len(patients)):
#     row = [pids[i]]
#     for j in range(len(patients)):
#         if i == j:
#             score = 1.0
#         else:
#             score = sims[(min(i, j), max(i, j))]
#         row.append('{:.4f}'.format(score))
#     print('\t'.join(row))
