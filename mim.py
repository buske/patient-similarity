#!/usr/bin/env python

"""

"""

__author__ = 'Orion Buske'

import os
import sys
import re
import logging

from collections import defaultdict

FREQUENCIES = {'very rare':  0.01, 
               'rare':       0.05, 
               'occasional': 0.075, 
               'frequent':   0.33, 
               'typical':    0.5, 
               'variable':   0.5, 
               'common':     0.75, 
               'hallmark':   0.9, 
               'obligate':   1.0}
fraction_frequency_re = re.compile(r'of|/')


class Disease:
    def __init__(self, db, id, name, phenotype_freqs):
        self.db = db
        self.id = id
        self.name = name
        self.phenotype_freqs = phenotype_freqs

class MIM:
    def __init__(self, filename):
        self.diseases = list(self.iter_diseases(filename))
        for i in range(5):
            logging.debug(self.diseases[i].__dict__)

    def __iter__(self):
        return iter(self.diseases)

    @classmethod
    def iter_disease_lines(cls, filename):
        with open(filename, encoding='utf-8') as ifp:
            cur_disease = None
            cur_lines = []
            for line in ifp:
                line = line.rstrip()
                tokens = line.split('\t')
                if len(tokens) == 1: continue
                disease = (tokens[0].strip(), tokens[1].strip())

                if disease == cur_disease:
                    cur_lines.append(tokens)
                else:
                    if cur_disease:
                        yield cur_disease, cur_lines

                    cur_lines = [tokens]
                    cur_disease = disease
            if cur_disease:
                yield cur_disease, cur_lines

    @classmethod
    def parse_frequency(cls, s, default=None):
        """Return float parsed frequency or default if problem or absent"""
        s = s.lower()
        if not s:
            freq = default
        elif s in FREQUENCIES:
            freq = FREQUENCIES[s]
        elif s.endswith('%'):
            s = s.replace('%', '')
            if '-' in s:
                # Average any frequency ranges
                low, high = s.split('-')
                freq = (float(low) + float(high)) / 2 / 100
            else:
                freq = float(s) / 100
        else:
            try:
                num, denom = fraction_frequency_re.split(s)
            except:
                logging.error("Error parsing frequency: {!r}".format(s))
                freq = default
            else:
                freq = float(num) / float(denom)

        return freq

    @classmethod
    def iter_diseases(cls, filename, default_freq=None):
        for disease, tokens_list in cls.iter_disease_lines(filename):
            db, id = disease
            raw_phenotypes = defaultdict(list)
            name = None
            for tokens in tokens_list:
                freq = cls.parse_frequency(tokens[8])
                hp_term = tokens[4].strip()
                raw_phenotypes[hp_term].append(freq)
                if not name:
                    name = tokens[2].strip()

            phenotype_freqs = {}
            for hp_term, freqs in raw_phenotypes.items():
                non_null = [x for x in freqs if x is not None]
                if non_null:
                    freq = sum(non_null) / len(non_null)
                else:
                    freq = default_freq

                phenotype_freqs[hp_term] = freq
                
            disease = Disease(db, id, name, phenotype_freqs)
            yield disease


def script(phenotype_filename, disease_gene_filename, out_hpo, out_genes, **kwargs):
    diseases = MIM(phenotype_filename)

    for disease in diseases:
        phenotypes = [('' if freq is None else '{:.6f}'.format(freq), hp) 
                      for (hp, freq) in disease.phenotype_freqs.items()]
        phenotypes.sort(reverse=True)
        for freq, hp in phenotypes:
            print('\t'.join(map(str, [disease.db, disease.id, hp, freq])))

    sys.exit(0)

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
    args = parse_args(args)

    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
