#!/usr/bin/env python

"""
Given all the information for a pair of patients, computer a match score, with gene and variant info
"""

import sys
import os
import logging
import csv

from collections import defaultdict

KO_THRESHOLD = 0.87  # between nonframeshift and splicing

def read_exomizer_vcf(filename):
    gene_scores = defaultdict(lambda: [None, [], None])  # gene -> (pheno, variant_scores, combined_score)
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
            combined = float(info['COMBINED_SCORE'])
            gene_scores[gene][0] = pheno
            gene_scores[gene][1].extend([geno] * n_alleles)
            gene_scores[gene][2] = combined

    for gene in gene_scores:
        # Sort and add a zero so every one has at least 2 values
        gene_scores[gene][1].append(0)
        gene_scores[gene][1].sort(reverse=True)

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
    sim_scores = defaultdict(list)  # id1 -> [(score, id2), ...]
    with open(filename) as ifp:
        for line in ifp:
            if not line or line.startswith('#'): continue
            tokens = line.strip().split('\t')
            p1, p2 = tokens[:2]
            p1 = ids.get(p1, p1)
            p2 = ids.get(p2, p2)
            score = float(tokens[2])
            sim_scores[p1].append((score, p2))
            sim_scores[p2].append((score, p1))

    logging.info('Read similarity scores for {} pairs'.format(len(sim_scores)))
    return sim_scores

def read_pheno_to_geno_file(filename):
    pheno_to_geno = {}
    with open(filename) as ifp:
        reader = csv.DictReader(ifp, delimiter=',')
        for row in reader:
            pid = row['Report ID'].strip()
            extern_id = row['Identifier'].strip()
            if not extern_id: continue
            assert extern_id not in pheno_to_geno
            pheno_to_geno[extern_id] = pid

    return pheno_to_geno

def score_gene(gene, p1, p2, patient_damages, sim_scores, 
               inheritance=None, control_damage=None):
    p1_damage = patient_damages[p1][gene]
    p2_damage = patient_damages[p2][gene]

    pheno_score = min(p1_damage[0], p2_damage[0])
    pheno_score = min((pheno_score + 0.1) / 0.7, 1)
    geno_1st = min(p1_damage[1][0], p2_damage[1][0])
    geno_2nd = min(p1_damage[1][1], p2_damage[1][1])

    score_dom = pheno_score * geno_1st
    score_rec = pheno_score * geno_2nd
    for other in set(patient_damages) - set([p1, p2]):
        other_damage = patient_damages[other].get(gene)
        if other_damage:
            pheno_similarity = max(sim_scores[p1, other], sim_scores[p2, other]) + 0.0001
            if other_damage[1][0] + 0.02 >= geno_1st:
                # hurt chances of dominant model
                score_dom *= pheno_similarity
            if other_damage[1][1] + 0.02 >= geno_2nd:
                # hurt chances of recessive model
                score_rec *= pheno_similarity

    if control_damage:
        # Use 1000gp as control damage              KOhom KOhet DMGhom DMGhet
        control_gene_damage = control_damage.get(gene, [0, 0, 0, 0])
        score_dom /= control_gene_damage[1] + 1
        if geno_1st < KO_THRESHOLD:
            score_dom /= control_gene_damage[3] + 1

        score_rec /= control_gene_damage[0] + 1
        if geno_2nd < KO_THRESHOLD:
            score_rec /= control_gene_damage[2] +1

    if not inheritance:
        return max(score_dom, score_rec)
    elif inheritance == 'AD':
        return score_dom
    elif inheritance == 'AR':
        return score_rec
    else:
        raise NotImplementedError('Unexpected inheritance: {}'.format(inheritance))

def average_score(gene, p1, p2, patient_damages, sim_scores,
                  inheritance=None, control_damage=None):
    p1_damage = patient_damages[p1][gene]
    p2_damage = patient_damages[p2][gene]
    return (p1_damage[2] + p2_damage[2]) / 2

def top_genes(p1, p2, patient_damages, sim_scores, inheritance=None,
              control_damage=None, method=None):
    p1_genes = patient_damages[p1]
    p2_genes = patient_damages[p2]
    shared_genes = set(p1_genes) & set(p2_genes)
    scores = []
    
    if method == 'avg':
        gene_scorer = average_score
    elif method == 'pc':
        gene_scorer = score_gene
    else:
        raise NotImplementedError('Unknown method: {}'.format(method))

    for gene in shared_genes:
        score = gene_scorer(gene, p1, p2, patient_damages, sim_scores,
                            inheritance=inheritance, control_damage=control_damage)
        scores.append((score, gene))

    #scores.sort(reverse=True)
    if scores:
        return max(scores)

def print_match(p1, p2, score, top):
    print('%s <-> %s: %.4f' % (p1, p2, score))
    for (score, gene) in top:
        score, inh = score
        if score >= 0.0001:
            print('    %.4f: %s (%s)' % (score, gene, inh))

def script(pheno_sim, exomiser_dir, inheritance=None, id_file=None,
           control_damage_file=None, method=None):
    pheno_scores = read_sim(pheno_sim)
    pair_scores = {}
    for p1, matches in pheno_scores.items():
        for score, p2 in matches:
            pair_scores[(p1, p2)] = score

    if control_damage_file:
        control_damage = read_gene_damages(control_damage_file)
        logging.info('Read control damage info: {}'.format(control_damage_file))
    else:
        control_damage = None

    if id_file:
        pheno_to_geno = read_pheno_to_geno_file(id_file)
        logging.info('Read patient IDs: {}'.format(id_file))
    else:
        pheno_to_geno = {}

    patient_damages = defaultdict(dict)
    for pid in pheno_scores:
        geno_id = pheno_to_geno.get(pid, pid)
        ezr_filename = os.path.join(exomiser_dir, geno_id + '.ezr')
        if os.path.isfile(ezr_filename):
            patient_damage = read_exomizer_vcf(ezr_filename)
            patient_damages[pid] = patient_damage
        else:
            logging.error('Missing EZR for: {}'.format(pid))

    logging.info('Read gene damage info for {} patients'.format(len(patient_damages)))
    logging.info('Using inheritance: {}'.format(inheritance))

    for p1 in sorted(pheno_scores):
        matches = pheno_scores[p1]
        #matches.sort(reverse=True)
        #for score, p2 in matches[:1]:
        score, p2 = max(matches)
        top = top_genes(p1, p2, patient_damages, pair_scores,
                        inheritance=inheritance, control_damage=control_damage,
                        method=method)
        top_score, top_gene = top if top else (float('nan'), '')
        print('{}\t{}\t{}\t{:.8f}'.format(p1, p2, top_gene, top_score))


def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()
    
    parser = ArgumentParser(description=description)
    parser.add_argument("pheno_sim")
    parser.add_argument("exomiser_dir")
    parser.add_argument('-I', '--inheritance', default=None,
                        choices=['AD', 'AR'])
    parser.add_argument('--method', default=None,
                        choices=['avg', 'pc'])
    parser.add_argument("--control-damage-file")
    parser.add_argument("--id-file", default=None)

    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()
        


    
