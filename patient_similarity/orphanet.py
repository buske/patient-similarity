#!/usr/bin/env python3

"""

"""

__author__ = 'Orion Buske'

import os
import sys
import logging

import xml.etree.ElementTree as ET

logger = logging.getLogger(__name__)

class Orphanet:
    def __init__(self, prevalence_filename, lookup_filename=None):
        lookup = None
        if lookup_filename:
            lookup = self.parse_lookup(lookup_filename)

        self.prevalence = self.parse_prevalence(prevalence_filename, lookup=lookup)

    def average_frequency(self):
        return sum(self.prevalence.values()) / len(self.prevalence)

    @classmethod
    def parse_lookup(cls, filename):
        tree = ET.parse(filename)
        root = tree.getroot()
        lookup = {}  # orphanet -> omim
        for disorder in root.findall('.//Disorder'):
            orphanum = disorder.find('OrphaNumber').text
            omim = None
            for ref in disorder.findall('.//ExternalReferenceList/ExternalReference'):
                source = ref.find('Source').text
                if source == 'OMIM':
                    omim = ref.find('Reference').text
                    break

            assert orphanum not in lookup
            if omim is not None:
                lookup[orphanum] = omim

        logger.info('Found {:d} Orphanet->OMIM entries'.format(len(lookup)))
        return lookup

    @classmethod
    def parse_prevalence(cls, filename, lookup=None):
        tree = ET.parse(filename)
        root = tree.getroot()
        prevalence = {}  # id -> prevalence

        prevalence_ids = {
            '12330': 2.5 / 1000,
            '12336': 7.5 / 10000,
            '12342': 2.5 / 10000,
            '12348': 5   / 100000,
            '12354': 5   / 1000000,
            '12360': 0.5 / 1000000,
            '12372': None,
            '12366': None
        }

        n_total = 0
        for disorder in root.findall('.//Disorder'):
            orphanum = disorder.find('OrphaNumber').text
            n_total += 1
            if lookup:
                try:
                    id = lookup[orphanum]
                except KeyError:
                    continue
            else:
                id = orphanum

            prevcls = disorder.find('ClassOfPrevalence')
            if not prevcls:
                continue

            previd = prevcls.get('id')
            prev = prevalence_ids[previd]
            if prev is None:
                continue

            prevalence[id] = prev

        logger.info('Found {:d} prevalences ({:d} dropped)'.format(len(prevalence), n_total - len(prevalence)))
        return prevalence


def script(lookup_filename, prevalence_filename):
    Orphanet(lookup_filename, prevalence_filename)

def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()

    parser = ArgumentParser(description=description)
    parser.add_argument('lookup_filename')
    parser.add_argument('prevalence_filename')
    return parser.parse_args()

def main(args=sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
