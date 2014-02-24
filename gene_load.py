#!/usr/bin/env python

"""

"""

__author__ = 'Orion Buske'

import os
import sys
import re
import math
import logging

from collections import defaultdict, Counter
from itertools import combinations, product

BROKEN_LOAD = 100
DAMAGED_LOAD = 10
NS_LOAD = 5  # should be based on background frequency

def pair_key(a, b):
    return (min(a, b), max(a, b))

def read_exomizer(filename):
    gene_scores = {}
    with open(filename) as ifp:
        for line in ifp:
            line = line.strip()
            if not line: continue
            tokens = line.split('\t')
            gene, pheno, geno, combined = tokens
            combined = float(combined)
            assert gene not in gene_scores
            gene_scores[gene] = combined
    return gene_scores

def compute_gene_weights(gene_scores):
    gene_weights = {}
    for gene, scores in gene_scores.items():
        if len(scores) >= 2:
            gene_weights[gene] = sum(scores)

    return gene_weights


def load_vcf_genes(filename):
    mut_re = re.compile(r'MUT=([^; \t,]+)')
    genes = Counter()  # Number of harmful mutation alleles observed so far
    n_dropped = 0
    with open(filename) as ifp:
        for line in ifp:
            line = line.strip()
            if line.startswith('#'): continue
            tokens = line.split('\t')
            info = tokens[7]
            gt = tokens[9].split(':', 1)[0]
            if not (gt == '0/1' or gt == '1/1'):
                n_dropped += 1
                continue

            mut_load = 1 if gt == '0/1' else 2
            match = mut_re.search(info).group(1)
            gene, mutation = match.split(':', 1)

            if 'wholegene' in mutation or \
                    ('p.' in mutation and (mutation.endswith('fs') or mutation.endswith('X'))):
                harmfulness = BROKEN_LOAD
            elif 'p.' in mutation and ('ins' in mutation or 'del' in mutation):
                harmfulness = DAMAGED_LOAD
            else:
                harmfulness = NS_LOAD

            genes[gene] += mut_load * harmfulness

    if n_dropped:
        logging.warning('Dropped {} variant lines with errors'.format(n_dropped))
    return genes
    

def read_load_data(filename):
    data = {}
    with open(filename) as ifp:
        for line in ifp:
            if line.startswith('#'): continue
            tokens = line.strip().split('\t')
            gene = tokens[0]
            data[gene] = dict(zip(['KK', 'K', 'DD', 'D', 'BB', 'B'], map(int, tokens[1:])))

    return data

def iter_vcf_entries(filename):
    with open(filename) as ifp:
        for line in ifp:
            if line.startswith('#'): continue
            id, vcf = line.strip().split('\t')
            assert ',' not in vcf
            yield vcf

def script(vcf_lookup, gene_load_data, exomizer_files, matrix=False):
    gene_load_data = read_load_data(gene_load_data)

    logging.info("Reading in genes from VCF files...")
    patient_gene_loads = {}  # id -> Counter
    for vcf in iter_vcf_entries(vcf_lookup):
        patient = os.path.basename(vcf).split('.')[0]
        gene_loads = load_vcf_genes(vcf)
        patient_gene_loads[patient] = gene_loads
    
    all_gene_scores = {}  # patient -> gene -> score
    patients = []
    gene_scores = defaultdict(list)  # gene -> scores
    for filename in exomizer_files:
        patient = os.path.basename(filename).split('.')[0]
        patient_scores = read_exomizer(filename)
        all_gene_scores[patient] = patient_scores
        patients.append(patient)
        for gene, score in patient_scores.items():
            gene_scores[gene].append(score)

    gene_weights = compute_gene_weights(gene_scores)

    pair_scores = {}  # pair -> (score, gene)
    for p1, p2 in combinations(patients, 2):
        s1 = all_gene_scores[p1]
        s2 = all_gene_scores[p2]
        shared_genes = set(s1) & set(s2)
        if shared_genes:
            scores = []
            for gene in shared_genes:
                gene_load = gene_load_data.get(gene)
                gene_weight = gene_weights[gene]
                load1 = patient_gene_loads[p1][gene]
                load2 = patient_gene_loads[p2][gene]
                min_load = min(load1, load2)
                if gene_load and (gene_load['KK'] or
                    (gene_load['K'] and min_load <= BROKEN_LOAD) or
                    (gene_load['DD'] and min_load < BROKEN_LOAD) or
                    (gene_load['D'] and min_load < DAMAGED_LOAD)):
                    continue

                pair_score = min(s1[gene], s2[gene])
                n_better = sum([1 for score in gene_scores[gene] if score >= pair_score]) - 1
                scores.append((pair_score / (n_better ** 2), gene))

            if scores:
                score, gene = max(scores)
                pair_scores[pair_key(p1, p2)] = (score, gene)

    if matrix:
        print('\t'.join(patients))
        for p1 in patients:
            p2_scores = []
            for p2 in patients:
                if p1 == p2:
                    score = 1
                else:
                    score, gene = pair_scores.get(pair_key(p1, p2), (0, ''))
                p2_scores.append('{:.4f}'.format(score))

            print('\t'.join([p1] + p2_scores))
    else:
        for p1, p2 in combinations(patients, 2):
            score, gene = pair_scores.get(pair_key(p1, p2), (0, ''))
            if score > 0:
                print('\t'.join([p1, p2, '{:.4f}'.format(score), gene]))
        
                
        
def parse_args():
    from argparse import ArgumentParser
    description = __doc__.strip()
    
    parser = ArgumentParser(description=description)
    parser.add_argument('vcf_lookup', metavar='VCF_LOOKUP')
    parser.add_argument('gene_load_data', metavar="GENE_LOADS",
                        help="1000gp/ESP6500 gene load data")
    parser.add_argument('exomizer_files', metavar='EXOMIZER_TSV', nargs='+')
    parser.add_argument('--matrix', action='store_true', default=False)
    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = parse_args()
    script(**vars(args))


    
