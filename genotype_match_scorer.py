#!/usr/bin/env python3

"""
Given patient similarity scores, incorporate genetic information and report candidate genes.
"""

import sys
import os
import logging
import csv

from math import log10
from collections import defaultdict
from bisect import bisect_left, bisect_right
from numpy import array
from scipy.stats import kendalltau

KO_THRESHOLD = 0.87  # between nonframeshift and splicing
EZR_SUFFIX = 'ezr'
EZR_FORMAT = 'ezr2'
N_MATCHES = 5  # phenotype matches/patient
N_GENES = 5  # genes/match
NAN = -1000

EZR_FORMATS = {
    'ezr1': {
        'gene': 'GENE',
        'pheno': 'PHENO_SCORE',
        'geno': 'VARIANT_SCORE',
        'combined': 'COMBINED_SCORE',
        'cadd': 'CADD',
        'cadd_phred': 'CADD_PHRED',
        },
    'ezr2': {
        'gene': 'EXOMISER_GENE',
        'pheno': 'EXOMISER_GENE_PHENO_SCORE',
        'geno': 'EXOMISER_VARIANT_SCORE',
        'combined': 'EXOMISER_GENE_COMBINED_SCORE',
        'cadd': 'CADD',
        'cadd_phred': 'CADD_PHRED',
        },
    'test': {
        'gene': 'GENE',
        'pheno': 'GENE_PHENO_SCORE',
        'geno': 'VARIANT_SCORE',
        'combined': 'GENE_COMBINED_SCORE',
        'cadd': 'CADD',
        'cadd_phred': 'CADD_PHRED',
        }
    }

def read_exomizer_vcf(filename, format=EZR_FORMAT):
    gene_scores = defaultdict(lambda: {'pheno': None, 'geno': [], 'combined': None, 'cadd': [], 'cadd_phred': []})
    col_names = EZR_FORMATS[format]
    with open(filename) as ifp:
        for line in ifp:
            if line.startswith('#'): continue
            try:
                tokens = line.rstrip('\n').split('\t')
                qual = float(tokens[5])
                if qual < 30: continue

                info = dict([part.split('=') for part in tokens[7].split(';') if '=' in part])
                gt = tokens[9]
                if gt == '-':
                    n_alleles = 1  # edge case for error with triallelics and exomizer
                else:
                    n_alleles = gt.count('1')
                if not 1 <= n_alleles <= 2:
                    logging.error('Unexpected gt: {!r}'.format(gt))
                    continue

                if col_names['pheno'] not in info:
                    col_names = EZR_FORMATS['test']
                    assert col_names['pheno'] in info
                    logging.error('SWITCHING TO TEST EZR FORMAT FOR FILE: {}!!!!'.format(filename))

                gene = info[col_names['gene']].upper()
                pheno = float(info[col_names['pheno']])
                geno = float(info[col_names['geno']])
                combined = float(info[col_names['combined']])
                cadd = float(info.get(col_names['cadd'], NAN))
                cadd_phred = float(info.get(col_names['cadd_phred'], 0))
                gene_scores[gene]['pheno'] = pheno
                gene_scores[gene]['geno'].extend([geno] * n_alleles)
                gene_scores[gene]['combined'] = combined
                gene_scores[gene]['cadd'].extend([cadd] * n_alleles)
                gene_scores[gene]['cadd_phred'].extend([cadd_phred] * n_alleles)
            except:
                logging.error('Error parsing line: {}'.format(line))
                raise

    for gene in gene_scores:
        # Sort and add a zero so every one has at least 2 values
        gene_scores[gene]['geno'].append(0)
        gene_scores[gene]['geno'].sort(reverse=True)
        gene_scores[gene]['cadd'].append(NAN)
        gene_scores[gene]['cadd'].sort(reverse=True)
        gene_scores[gene]['cadd_phred'].append(0)
        gene_scores[gene]['cadd_phred'].sort(reverse=True)

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
        header = ifp.readline()
        assert header.startswith('A\tB\t')
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

def scale_pheno_score(pheno_scale, score):
    if pheno_scale == 'ezr1':
        return min((score + 0.1) / 1.1, 1)
    elif pheno_scale == 'ezr2':
        return score
    else:
        raise NotImplementedError('Unknown pheno scale: {}'.format(pheno_scale))

def pc_score(gene, p1, p2, patient_damages, sim_scores, pheno_scale=None,
              inheritance=None, control_damage=None, *args, **kwargs):
    p1_damage = patient_damages[p1][gene]
    p2_damage = patient_damages[p2][gene]

    pheno_score = min(p1_damage['pheno'], p2_damage['pheno'])
    pheno_score = scale_pheno_score(pheno_scale, pheno_score)

    geno_1st = min(p1_damage['geno'][0], p2_damage['geno'][0])
    geno_2nd = min(p1_damage['geno'][1], p2_damage['geno'][1])

    score_dom = pheno_score * geno_1st
    score_rec = pheno_score * geno_2nd
    for other in set(patient_damages) - set([p1, p2]):
        other_damage = patient_damages[other].get(gene)
        if other_damage:
            pheno_similarity = max(sim_scores[p1, other], sim_scores[p2, other]) + 0.001
            if other_damage['geno'][0] + 0.01 >= geno_1st:
                # hurt chances of dominant model
                score_dom *= pheno_similarity
            if other_damage['geno'][1] + 0.01 >= geno_2nd:
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

def pc_score_test(gene, p1, p2, patient_damages, sim_scores, pheno_scale=None,
              inheritance=None, control_damage=None, *args, **kwargs):
    p1_damage = patient_damages[p1][gene]
    p2_damage = patient_damages[p2][gene]

    pheno_score = min(p1_damage['pheno'], p2_damage['pheno'])
    pheno_score = scale_pheno_score(pheno_scale, pheno_score)

    geno_1st = min(p1_damage['geno'][0], p2_damage['geno'][0])
    geno_2nd = min(p1_damage['geno'][1], p2_damage['geno'][1])

    score_dom = pheno_score * geno_1st
    score_rec = pheno_score * geno_2nd
    for other in set(patient_damages) - set([p1, p2]):
        other_damage = patient_damages[other].get(gene)
        if other_damage:
            pheno_similarity = min(max(sim_scores[p1, other], sim_scores[p2, other]) + 0.01, 1)
            if other_damage['geno'][0] >= geno_1st:
                # hurt chances of dominant model
                score_dom *= pheno_similarity
            if other_damage['geno'][1] >= geno_2nd:
                # hurt chances of recessive model
                score_rec *= pheno_similarity

    if control_damage:
        # Use 1000gp as control damage              KOhom KOhet DMGhom DMGhet
        control_gene_damage = control_damage.get(gene, [0, 0, 0, 0])
        score_dom *= 0.5 ** (control_gene_damage[1])
        if geno_1st < KO_THRESHOLD:
            score_dom *= 0.5 ** (control_gene_damage[3])

        score_rec *= 0.5 ** (control_gene_damage[0])
        if geno_2nd < KO_THRESHOLD:
            score_rec *= 0.5 ** (control_gene_damage[2])

    if not inheritance:
        return max(score_dom, score_rec)
    elif inheritance == 'AD':
        return score_dom
    elif inheritance == 'AR':
        return score_rec
    else:
        raise NotImplementedError('Unexpected inheritance: {}'.format(inheritance))

def pc_cadd_score(gene, p1, p2, patient_damages, sim_scores, cadd_distributions=None, control_damage=None,
               inheritance=None, pheno_scale=None, *args, **kwargs):
    p1_damage = patient_damages[p1][gene]
    p2_damage = patient_damages[p2][gene]

    pheno_score = min(p1_damage['pheno'], p2_damage['pheno'])
    pheno_score = scale_pheno_score(pheno_scale, pheno_score)

    geno = [0, 0]
    inh_scores = [0, 0]
    for i in [0, 1]:
        geno[i] = min(p1_damage['cadd'][i], p2_damage['cadd'][i])
        cadd_phred = min(p1_damage['cadd_phred'][i], p2_damage['cadd_phred'][i])
        cadd_p = 10 ** (-cadd_phred / 10)
        assert 0 <= cadd_p <= 1
        inh_scores[i] = pheno_score * (1 - cadd_p)

    for other in set(patient_damages) - set([p1, p2]):
        other_damage = patient_damages[other].get(gene)
        if other_damage:
            pheno_similarity = max(sim_scores[p1, other], sim_scores[p2, other]) + 0.001
            for i in [0, 1]:
                if other_damage['cadd'][i] + 0.01 >= geno[i]:
                    # hurt chances of inheritance model
                    inh_scores[i] *= pheno_similarity

    if control_damage:
        raise NotImplementedError('control_damage')

    if cadd_distributions:
        raise NotImplementedError('cadd_distributions')

    if not inheritance:
        return max(inh_scores)
    elif inheritance == 'AD':
        return inh_scores[0]
    elif inheritance == 'AR':
        return inh_scores[1]
    else:
        raise NotImplementedError('Unexpected inheritance: {}'.format(inheritance))



def average_score(gene, p1, p2, patient_damages, *args, **kwargs):
    p1_damage = patient_damages[p1][gene]
    p2_damage = patient_damages[p2][gene]
    return (p1_damage['combined'] + p2_damage['combined']) / 2


def cadd_score(gene, p1, p2, patient_damages, sim_scores, cadd_distributions,
               inheritance=None, pheno_scale=None, *args, **kwargs):
    p1_damage = patient_damages[p1][gene]
    p2_damage = patient_damages[p2][gene]

    pheno_score = min(p1_damage['pheno'], p2_damage['pheno'])
    pheno_score = scale_pheno_score(pheno_scale, pheno_score)

    geno = [0, 0]
    inh_scores = [0, 0]
    for i in [0, 1]:
        geno[i] = min(p1_damage['cadd'][i], p2_damage['cadd'][i])
        cadd_phred = min(p1_damage['cadd_phred'][i], p2_damage['cadd_phred'][i])
        cadd_p = 10 ** (-cadd_phred / 10)
        assert 0 <= cadd_p <= 1
        inh_scores[i] = pheno_score * (1 - cadd_p)

    for other in set(patient_damages) - set([p1, p2]):
        other_damage = patient_damages[other].get(gene)
        if other_damage:
            pheno_similarity = min(max(sim_scores[p1, other], sim_scores[p2, other]) + 0.5, 1)
            for i in [0, 1]:
                if other_damage['cadd'][i] >= geno[i]:
                    # hurt chances of inheritance model
                    inh_scores[i] *= pheno_similarity

    if cadd_distributions:
        for i in [0, 1]:
            if gene in cadd_distributions[i]:
                distribution = cadd_distributions[i][gene]
                n_worse = len(distribution) - bisect_left(distribution, geno[i])
                inh_scores[i] *= 0.5 ** n_worse

    if not inheritance:
        return max(inh_scores)
    elif inheritance == 'AD':
        return inh_scores[0]
    elif inheritance == 'AR':
        return inh_scores[1]
    else:
        raise NotImplementedError('Unexpected inheritance: {}'.format(inheritance))


def get_scored_genes(p1, p2, patient_damages, sim_scores, inheritance=None,
                     control_damage=None, method=None, cadd_distributions=None,
                     pheno_scale=None):
    p1_genes = patient_damages[p1]
    p2_genes = patient_damages[p2]
    shared_genes = set(p1_genes) & set(p2_genes)
    scores = []
    
    if method == 'avg':
        gene_scorer = average_score
    elif method == 'pc':
        gene_scorer = pc_score
    elif method == 'pc-cadd':
        gene_scorer = pc_cadd_score
    elif method == 'pc-test':
        gene_scorer = pc_score_test
    elif method == 'cadd':
        gene_scorer = cadd_score
    else:
        raise NotImplementedError('Unknown method: {}'.format(method))

    for gene in shared_genes:
        score = gene_scorer(gene, p1, p2, patient_damages, sim_scores=sim_scores,
                            inheritance=inheritance, control_damage=control_damage,
                            cadd_distributions=cadd_distributions,
                            pheno_scale=pheno_scale)
        scores.append((score, gene))

    return scores

def read_solution_genes(filename):
    solutions = {}
    with open(filename) as ifp:
        for line in ifp:
            patient, genes = line.rstrip('\n').split('\t')
            genes = genes.strip().upper().split(',')
            assert patient not in solutions
            if genes:
                solutions[patient] = set(genes)

    return solutions

def load_cadd_distribution(filename):
    distribution = {}  # gene -> distribution
    with open(filename) as ifp:
        for line in ifp:
            tokens = line.rstrip('\n').replace('nan', str(NAN)).split('\t')
            gene = tokens[0]
            gene_distribution = sorted(map(float, tokens[1:]))
            distribution[gene] = gene_distribution

    logging.info('Loaded CADD distributions for {} genes, e.g.:\n{}\t{}...'.format(len(distribution), gene, distribution[gene][1:10]))

    return distribution

def load_cadd_distributions(base):
    return [load_cadd_distribution(base + '.0.txt'),
            load_cadd_distribution(base + '.1.txt')]

def script(pheno_sim, exomiser_dir, inheritance=None, id_file=None,
           control_damage_file=None, solution_gene_file=None, pheno_scale=None,
           cadd_base=None, method=None, ezr_suffix=EZR_SUFFIX, ezr_format=EZR_FORMAT):
    logging.info('Using inheritance: {}'.format(inheritance))
    logging.info('Using gene_pheno_score scaling {} such that 0.5 -> {}'.format(pheno_scale, scale_pheno_score(pheno_scale, 0.5)))

    pheno_scores = read_sim(pheno_sim)
    pair_scores = {}
    for p1, matches in pheno_scores.items():
        for score, p2 in matches:
            pair_scores[(p1, p2)] = score

    if cadd_base:
        cadd_distributions = load_cadd_distributions(cadd_base)
    else:
        cadd_distributions = None

    if control_damage_file:
        control_damage = read_gene_damages(control_damage_file)
        logging.info('Read control damage info: {}'.format(control_damage_file))
    else:
        control_damage = None

    if solution_gene_file:
        solution_genes = read_solution_genes(solution_gene_file)
        n_overlap = len(set(pheno_scores).intersection(solution_genes))
        logging.info('Read candidate genes for {} patients ({} with phenotypes)'.format(len(solution_genes), n_overlap))
    else:
        solution_genes = None

    if id_file:
        pheno_to_geno = read_pheno_to_geno_file(id_file)
        logging.info('Read patient IDs: {}'.format(id_file))
    else:
        pheno_to_geno = {}

    patient_damages = defaultdict(dict)
    incomplete_patients = set()
    for pid in pheno_scores:
        geno_id = pheno_to_geno.get(pid, pid)
        ezr_filename = os.path.join(exomiser_dir, geno_id + '.' + ezr_suffix)
        if os.path.isfile(ezr_filename):
            try:
                patient_damage = read_exomizer_vcf(ezr_filename, format=ezr_format)
            except:
                logging.error('Encountered error reading file: {}'.format(ezr_filename))
                raise

            patient_damages[pid] = patient_damage
        else:
            logging.error('Missing EZR for: {}'.format(pid))
            incomplete_patients.add(pid)


    logging.info('Read gene damage info for {} patients with phenotypes'.format(len(patient_damages)))

    for p1 in sorted(pheno_scores):
        if p1 in incomplete_patients or p1 not in solution_genes: continue

        matches = pheno_scores[p1]
        matches.sort(reverse=True)
        match_i = 0
        for pheno_score, p2 in matches:
            if p2 in incomplete_patients: continue
            if match_i >= N_MATCHES: break
            match_i += 1

            scored_genes = get_scored_genes(p1, p2, patient_damages, pair_scores, 
                                            inheritance=inheritance, control_damage=control_damage,
                                            method=method, cadd_distributions=cadd_distributions,
                                            pheno_scale=pheno_scale)
            scored_genes.sort(reverse=True)
            top_n = scored_genes[:N_GENES]
            output = [p1, p2, pheno_score]
            if solution_genes:
                # Report rank and score of top true/candidate gene
                true_genes = solution_genes[p1]
                hit = [((i + 1), gene) for (i, (_, gene)) in enumerate(scored_genes) if gene in true_genes]
                hit_rank, hit_gene = hit[0] if hit else ('NA', 'NA')

                output.append('{}:{}'.format(hit_gene, hit_rank))

            # Report scores of top genes
            if top_n:
                for top_score, top_gene in top_n:
                    output.append('{}:{}'.format(top_gene, top_score))
            else:
                output.append('NA')

            print('\t'.join(map(str, output)))


def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()
    
    parser = ArgumentParser(description=description)
    parser.add_argument("pheno_sim")
    parser.add_argument("exomiser_dir")
    parser.add_argument('-I', '--inheritance', default=None,
                        choices=['AD', 'AR'])
    parser.add_argument('--method', default=None,
                        choices=['avg', 'pc', 'cadd', 'pc-cadd', 'pc-test'])
    parser.add_argument("--control-damage-file")
    parser.add_argument("--cadd-base")
    parser.add_argument("--solution-gene-file")
    parser.add_argument("--id-file", default=None)
    parser.add_argument("--ezr-suffix", default=EZR_SUFFIX)
    parser.add_argument("--ezr-format", default=EZR_FORMAT)
    parser.add_argument("--pheno-scale", default=None,
                        choices=['ezr1', 'ezr2'])
    parser.add_argument("--log", default='INFO')

    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = vars(parse_args(args))
    logging.basicConfig(level=args.pop('log'))
    script(**args)

if __name__ == '__main__':
    sys.exit(main())
        


    
