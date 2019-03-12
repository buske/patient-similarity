#!/usr/bin/env python3

"""
Helper script to generate disease-prototypes of common HPO terms and genes
"""

__author__ = 'Orion Buske'

import os
import sys
import re
import logging

from collections import defaultdict

from patient_similarity import Diseases

logger = logging.getLogger(__name__)


def script(phenotype_filename, disease_gene_filename, out_hpo, out_genes, **kwargs):
    diseases = Diseases(phenotype_filename, db='ORPHANET')

    disease_genes = defaultdict(set)
    with open(disease_gene_filename) as ifp:
        for line in ifp:
            line = line.strip()
            if not line or line.startswith('#'): continue
            tokens = line.split('\t')
            if len(tokens) == 1: continue
            disease_id, gene_id, gene_name = tokens
            disease_genes[disease_id].add(gene_name)

    # Generate prototypical patients
    with open(out_hpo, 'w') as out_hpo, open(out_genes, 'w') as out_genes:
        for disease in diseases:
            hp_terms = []
            for hp_term, freq in disease.phenotype_freqs.items():
                if freq is None or freq >= 0.25:
                    hp_terms.append(hp_term)

            disease_id = '{}:{}'.format(disease.db, disease.id)
            genes = disease_genes[disease_id]
            if len(hp_terms) >= 5 and genes:
                disease_id = 'PROTO:{}'.format(disease_id)
                print('{}\t{}'.format(disease_id, ';'.join(hp_terms)), file=out_hpo)
                print('{}\t{}'.format(disease_id, ';'.join(genes)), file=out_genes)


def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()

    parser = ArgumentParser(description=description)
    parser.add_argument('phenotype_filename', metavar='phenotype_annotations.tab')
    parser.add_argument('disease_gene_filename', metavar='diseases_to_genes.txt')
    parser.add_argument('out_hpo', metavar='out.hpo')
    parser.add_argument('out_genes', metavar='out.genes')

    return parser.parse_args(args)


def main(args=sys.argv[1:]):
    logging.basicConfig(level='INFO')
    args = parse_args(args)
    script(**vars(args))


if __name__ == '__main__':
    sys.exit(main())
