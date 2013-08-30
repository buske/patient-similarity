#!/usr/bin/env python

from all_data import *
from phenotype import *
import numpy as np
import sim
import argparse
import sys
import pprint

parser = argparse.ArgumentParser(description='Reads in a list of patients. Each line of input should contain a list of HPO terms separated by spaces describing a single patient. The algorithm computes the similarity of the first patient to each of the remaining patients. It takes a bit of time to load up some pre-computed values.')
parser.add_argument('-x', dest='x', action='store_true', help='If this flag is on, the algorithm uses a slightly better bound on the conditional probabilities. See the readme for details.')
parser.add_argument('-n', dest='n', action='store', type=int, default = -1, help='Number of patients to process. If not provided, processes until an empty line is provided.')
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='Input file. If not provided, reads from stdin.')
args = parser.parse_args()

def load_patient(s):
	return Patient([code_to_term[x] for x in s.split()])

lines = args.infile
p1 = load_patient(lines.readline())

if args.n >= 0:
	for i in range(args.n - 1):
		print sim.sim(p1, load_patient(lines.readline()), args.x)
else:
	while True:
		l = lines.readline()
		if len(l) < 2:
			break
		print sim.sim(p1,load_patient(l), args.x)
