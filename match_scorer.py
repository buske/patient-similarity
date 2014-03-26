#!/usr/bin/env python

"""
Given all the information for a pair of patients, computer a match score, with gene and variant info
"""

import sys
import os
import math
import logging
from collections import defaultdict
from itertools import combinations, product
from numpy import array

args = sys.argv[1:]

KO_THRESHOLD = 0.87  # between nonframeshift and splicing

def read_ids(filename):
    pheno_ids = {}
    geno_ids = {}
    ids = set()
    with open(filename) as ifp:
        for line in ifp:
            line = line.strip()
            if not line or line.startswith('#'): continue
            id, pheno_id, geno_id = line.split('\t')[:3]
            if pheno_id == '-' or geno_id == '-':
                # Missing phenotype or genotype data
                continue

            assert pheno_id not in pheno_ids
            assert geno_id not in geno_ids
            assert id not in ids
            pheno_ids[pheno_id] = id
            geno_ids[geno_id] = id
            ids.add(id)
            

    logging.info('Found {} complete patients'.format(len(ids)))
    return pheno_ids, geno_ids, ids


def read_causal_genes(filename):
    causal_genes = {}  # -> (id, gene)
    with open(filename) as ifp:
        for line in ifp:
            line = line.strip()
            if not line or line.startswith('#'): continue
            tokens = line.split('\t')
            if len(tokens) > 1:
                id, gene = tokens
                causal_genes[id] = gene            

    logging.info('Found causal gene for {} patients'.format(len(causal_genes)))
    return causal_genes


def read_exomizer_vcf(filename):
    gene_scores = defaultdict(lambda: [None, []])  # gene -> pheno, scores
    with open(filename) as ifp:
        for line in ifp:
            line = line.strip()
            if not line or line.startswith('#'): continue

            tokens = line.split('\t')
            info = dict([part.split('=') for part in tokens[7].split(';')])
            gt = tokens[9]
            if gt == '-':
                n_alleles = 1  # edge case for error with triallelics and exomizer
            else:
                n_alleles = gt.count('1')
            if not 1 <= n_alleles <= 2:
                logging.error('Unexpected gt: {!r}'.format(gt))
                continue

            gene = info['GENE']
            pheno = float(info['PHENO_SCORE'])
            geno = float(info['VARIANT_SCORE'])

            gene_scores[gene][0] = pheno
            gene_scores[gene][1].extend([geno] * n_alleles)

    for gene in gene_scores:
        # Add a zero so every one has at least 2 values
        gene_scores[gene][1] = array(sorted(gene_scores[gene][1] + [0]))

    return dict(gene_scores)

def read_gene_damages(filename):
    gene_scores = {}
    with open(filename) as ifp:
        for line in ifp:
            line = line.strip()
            if not line or line.startswith('#'): continue
            tokens = line.split('\t')
            gene = tokens[0]
            damages = list(map(int, tokens[1:5]))
            assert gene not in gene_scores
            gene_scores[gene] = damages

    logging.info('Read gene damage stats for {} genes'.format(len(gene_scores)))
    return gene_scores

def read_sim(filename, ids={}):
    sim_scores = {}  # (id1, id2) -> score
    with open(filename) as ifp:
        for line in ifp:
            if not line or line.startswith('#'): continue
            tokens = line.strip().split('\t')
            p1, p2 = tokens[:2]
            p1 = ids.get(p1, p1)
            p2 = ids.get(p2, p2)
            scores = list(map(float, tokens[2:]))
            score, shared_ic = scores
            sim_scores[(p1, p2)] = score

    logging.info('Read similarity scores for {} pairs'.format(len(sim_scores)))
    return sim_scores

def score_gene(gene, p1, p2, control_damages, patient_damages, sim_scores):
    p1_damage = patient_damages[p1][gene]
    p2_damage = patient_damages[p2][gene]

    # Use 1000gp as control damage             KOhom KOhet DMGhom DMGhet
    control_damage = control_damages.get(gene, [0, 0, 0, 0])

    gene_pheno = min(p1_damage[0], p2_damage[0])
    gene_geno1 = min(p1_damage[1][-1], p2_damage[1][-1])
    gene_geno2 = min(p1_damage[1][-2], p2_damage[1][-2])

    score_dom = 1
    score_rec = 1
    for other in set(patient_damages) - set([p1, p2]):
        other_damage = patient_damages[other].get(gene)
        if other_damage:
            pheno_similarity = min(sim_scores[p1, other], sim_scores[p2, other]) + 0.001
            if other_damage[1][-1] >= gene_geno1:
                # hurt chances of dominant model
                score_dom *= pheno_similarity
            if other_damage[1][-2] >= gene_geno2:
                # hurt chances of recessive model
                score_rec *= pheno_similarity

    score_dom /= control_damage[1] + 1
    if gene_geno1 < KO_THRESHOLD:
        score_dom /= control_damage[3] + 1

    score_rec /= control_damage[0] + 1
    if gene_geno2 < KO_THRESHOLD:
        score_rec /= control_damage[2] +1

    return max((score_dom, 'D'), (score_rec, 'r'))


def top_genes(p1, p2, control_damages, patient_damages, sim_scores):
    p1_genes = patient_damages[p1]
    p2_genes = patient_damages[p2]
    shared_genes = set(p1_genes) & set(p2_genes)
    scores = []
    for gene in shared_genes:
        score = score_gene(gene, p1, p2, control_damages, patient_damages, sim_scores)
        scores.append((score, gene))

    scores.sort(reverse=True)
    return scores[:5]

def print_match(p1, p2, score, top):
    print('%s <-> %s: %.4f' % (p1, p2, score))
    for (score, gene) in top:
        score, inh = score
        if score >= 0.0001:
            print('    %.4f: %s (%s)' % (score, gene, inh))

def script(id_lookup, pheno_sim, causal_file, gene_damage_file, exomizer_files):
    pheno_ids, geno_ids, ids = read_ids(id_lookup)
    causal_genes = read_causal_genes(causal_file)
    control_damages = read_gene_damages(gene_damage_file)

    patient_damages = defaultdict(dict)
    for filename in exomizer_files:
        patient_damage = read_exomizer_vcf(filename)
        geno_id = os.path.basename(filename).split('.')[0]
        id = geno_ids[geno_id]
        patient_damages[id] = patient_damage

    logging.info('Read gene damage info for {} patients'.format(len(patient_damages)))
    scores = read_sim(pheno_sim, ids=pheno_ids)

    ids = ["UDP_1019", "UDP_2058", "380_120891B", "UDP_3384", "UDP_2803", "UDP_5730", "464_MT0003", "UDP_4306", "UDP_5316", "474_BC0006", "474_BC0005"]
    for p1 in ids:
        # p1 = '174_10-462'
        matches = [(scores.get((p1, p), 0), p) for p in ids]
        matches.sort(reverse=True)

        for score, p2 in matches[:2]:
            top = top_genes(p1, p2, control_damages, patient_damages, scores)
            print_match(p1, p2, score, top)
            # # print('\t'.join(map(str, [p1, p2, z1, z2, pair_sim[1], p1_ic, p2_ic, len(shared_genes)])))
            
            # gene_scores = [(score_gene(gene, p1, p2, control_damages, patient_damages), gene) for gene in shared_genes]
            # gene_scores.sort(reverse=True)
            # print('\t'.join(map(str, [p1, p2, min_z, '\t'.join(map(str, gene_scores[:5]))])))
            # For each shared mutated gene, compute score for dominant and recessive modes
            # depending on damage in the gene within two patients and all others
            

            #p1_damage[1][-2:]
            #p2_damage[1][-2:]

            # if p1_damage[0] > 2 or p2_damage[0] > 2 or control_damage[0] > 0: continue
            # if (p1_damage[0] == 0 and p1_damage[1] > 3) or (p2_damage[0] == 0 and p2_damage[1] > 3): continue
            # if (p1_damage[0] <= 1 and p2_damage[0] <= 1) and control_damage[1] > 0: continue
            # if (p1_damage[0] == 0 or p2_damage[0] == 0) and control_damage[2] > 0: continue
            
            #print('{}\t{}\t{}\t{}'.format(gene, p1_damage, p2_damage, case_damages))

        # break
    

def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()
    
    parser = ArgumentParser(description=description)
    parser.add_argument("id_lookup")
    parser.add_argument("pheno_sim")
    parser.add_argument("causal_file", metavar="causal_genes")
    parser.add_argument("gene_damage_file", metavar="gene_damages")
    parser.add_argument("exomizer_files", metavar="exomizer", nargs='+')

    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()
        


    
