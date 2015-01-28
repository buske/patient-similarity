#!/usr/bin/env python3

"""
Lookup information about a list of HPO terms, specified on stdin.
"""

__author__ = 'Orion Buske (buske@cs.toronto.edu)'

import sys
import logging

from hpo import HPO

def script(hpo_filename):
    hpo = HPO(hpo_filename)

    for line in sys.stdin:
        term = line.strip()
        if not term: continue

        assert term.startswith('HP:')
        print('[{0.id}] {0.name}'.format(hpo[term]))

def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()
    
    parser = ArgumentParser(description=description)
    parser.add_argument('hpo_filename', metavar='hp.obo')
    parser.add_argument('--log', dest='loglevel', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default='WARNING')

    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = parse_args(args)
    logging.basicConfig(level=args.__dict__.pop('loglevel'))
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
