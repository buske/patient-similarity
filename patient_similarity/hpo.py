#!/usr/bin/env python3

"""
Module for processing the Human Phenotype Ontology.

Provides a Class interface for interacting with the HPO:
> hpo = HPO('hpo.obo')
> hpo.root
HP:0000001

Restrict to a sub-branch of the HPO:
> hpo.filter_to_descendants('HP:0000118')
> hpo.root
HP:0000118

Lookup HPNode by id or alt_id:
> term = hpo['HP:0000801']
> term.name
'Nephrotic syndrome'
> term.id  # canonical id
'HP:0000100'

Can get parents or children HPNodes:
> term.parents
{HP:0012211}
> term.children
{HP:0012589, HP:0012588, HP:0008695, HP:0008677}
"""

__author__ = 'Orion Buske (buske@cs.toronto.edu)'

import sys

assert sys.version_info >= (3, 0), 'Python 3 is required'

import os
import re
import logging

logger = logging.getLogger(__name__)

class HPError(Exception):
    pass

class HPObsoleteError(HPError):
    pass

def get_descendants(root, acc=None):
    """Add to acc all descendants of root"""
    if acc is None:
        acc = set()
    acc.add(root)
    for child in root.children:
        get_descendants(child, acc)
    return acc

def get_ancestors(root, acc=None):
    if acc is None:
        acc = set()
    acc.add(root)
    for parent in root.parents:
        get_ancestors(parent, acc)
    return acc


class HPNode(object):
    """HPO graph node

    Attributes:
    id
    name
    parents (HP objects, filled in by HPO)
    children (HP objects, filled in by HPO)
    alts (HP terms)
    """
    def __init__(self, lines):
        self.parents = set()
        self.children = set()
        self.alts = set()
        self._parent_hps = set()
        for line in lines:
            line = line.strip()
            field, value = line.split(': ', 1)
            if field == 'is_obsolete':
                raise HPObsoleteError()
            elif field == 'id':
                assert value.startswith('HP:') and value[-1].isdigit(), value
                self.id = value
            elif field == 'name':
                self.name = value
            elif field == 'alt_id':
                assert value.startswith('HP:') and value[-1].isdigit(), value
                self.alts.add(value)
            elif field == 'is_a':
                hp = value.split('!')[0].strip()
                assert hp.startswith('HP:') and hp[-1].isdigit(), value
                self._parent_hps.add(hp)

        try:
            assert self.id
            assert self.name
            if self.id != 'HP:0000001':
                assert self._parent_hps
        except:
            logger.error("Error parsing TERM:\n{}".format('\n'.join(lines)))
            raise

    def __str__(self):
        return str(self.id)

    def __repr__(self):
        return str(self)

    def link(self, hps):
        """Link to objects for parents and children, given lookup dict"""
        for hp in self._parent_hps:
            parent = hps[hp]
            self.parents.add(parent)
            parent.children.add(self)

    def is_root(self):
        return len(self._parent_hps) == 0

    def ancestors(self):
        return get_ancestors(self)

    def descendants(self):
        return get_descendants(self)

def _iter_hp_terms(reader):
    term_lines = None
    for line in reader:
        line = line.strip()
        if not line: continue
        if line == '[Term]':
            if term_lines:
                yield term_lines

            term_lines = []
        else:
            if term_lines is not None:
                term_lines.append(line)

    if term_lines:
        yield term_lines


class HPO(object):
    """HPO graph

    Attributes:
    version
    hps: {hp -> HP}
    root: HP
    """
    def __init__(self, filename, new_root=None):
        """Load the HPO ontology specified in the OBO file 'filename'.

        If new_root is specified (an HPO term ID), the HPO will be truncated to
        only include that node and descendants.
        """
        self.hps = {}
        self.root = None

        logger.info("Parsing graph...")
        roots = []
        with open(filename, encoding='utf-8') as ifp:
            version_str = ifp.readline().strip()
            assert version_str.startswith('format-version')
            self.version = version_str.split(': ')[1]
            logger.info("HPO version: {}".format(self.version))

            for lines in _iter_hp_terms(ifp):
                try:
                    hp = HPNode(lines)
                except HPError:
                    continue

                if hp.is_root():
                    roots.append(hp)
                # Relate all alt HP ids to object
                self.hps[hp.id] = hp
                for hpid in hp.alts:
                    assert hpid not in self.hps
                    self.hps[hpid] = hp

        # Connect network of HP objects
        nodes = set(self.hps.values())
        for node in nodes:
            node.link(self.hps)

        logger.info("Found {:d} HP nodes ({:d} terms) in graph".format(len(nodes), len(self.hps)))

        logger.debug("Here are 5:")
        for i, k in zip(list(range(5)), nodes):
            logger.debug("  {:d}: {}".format(i, k))

        if len(roots) == 1:
            self.root = roots[0]
        else:
            logger.warning("Warning: found {:d} root nodes, leaving root as None".format(len(roots)))
            self.root = None

        if new_root:
            assert new_root in self.hps
            self.filter_to_descendants(new_root)


    def filter_to_descendants(self, root_hp):
        root = self.hps[root_hp]

        safe_nodes = get_descendants(root)
        logger.info("Filtering to the {:d} nodes descendant of {} ({})...".format(len(safe_nodes), root_hp, root.name))

        hps = {}
        safe_hps = set()
        for node in safe_nodes:
            safe_hps.add(node.id)
            safe_hps.update(node.alts)
            if node.id in hps:
                logger.warning('Found duplicate of HP term:' + node.id)

            hps[node.id] = node
            for alt_hp in node.alts:
                if alt_hp in hps:
                    logger.warning('Found duplicate of HP term:' + alt_hp)
                else:
                    hps[alt_hp] = node

        # Reset all connections in network
        for node in safe_nodes:
            node.parents.clear()
            node.children.clear()
            node._parent_hps.intersection_update(safe_hps)

        # Re-link
        for node in safe_nodes:
            node.link(self.hps)

        # Replace attributes
        self.root = root
        self.hps = hps

    def __getitem__(self, key):
        return self.hps[key]

    def __iter__(self):
        return iter(set(self.hps.values()))

    def __len__(self):
        return len(set(self.hps.values()))

    def descendant_terms(self, root_hp):
        root = self.hps[root_hp]

        descendants = get_descendants(root)
        terms = set()
        for node in descendants:
            terms.add(node.id)
            terms.update(node.alts)

        return terms


def script(hpo_filename):
    hpo = HPO(hpo_filename)
    hpo.filter_to_descendants('HP:0000118')

def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()

    parser = ArgumentParser(description=description)
    parser.add_argument('hpo_filename', metavar='hp.obo')
    return parser.parse_args()

def main(args=sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
