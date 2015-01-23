#!/usr/bin/env python

"""
Module for helper functions and other things that don't fit neatly elsewhere.
"""

import sys

import numpy as np
from numpy import linalg as LA

import pymatgen as mg
from pymatgen.symmetry.finder import SymmetryFinder

from splpy.resio import Res


def find_or_create(collection, query, update, **kwargs):
    return collection.find_and_modify(query, {"$setOnInsert": update}, upsert=True, new=True, **kwargs)


class OrderedPair(object):
    """
    Object for representing a pair of objects (usually numbers) that are
    accessed in order using the .first and .second attributes

    """

    def __init__(self, x0=None, x1=None):
        self._first = None
        self._second = None
        self.set(x0, x1)

    @property
    def first(self):
        return self._first

    @property
    def second(self):
        return self._second

    def set(self, x0, x1):
        if x0 < x1:
            self._first = x0
            self._second = x1
        else:
            self._first = x1
            self._second = x0

    def __str__(self):
        return "{}~{}".format(self.first, self.second)

    def __hash__(self):
        return hash(self.first) ^ hash(self.second)

    def __cmp__(self, other):
        """Compare."""
        if other is None or not isinstance(other, OrderedPair):
            return 1
        start_cmp = cmp(self.first, other.first)
        if 0 != start_cmp:
            return start_cmp
        else:
            return cmp(self.second, other.second)

    @classmethod
    def from_string(cls, s):
        values = s.split("~")
        if len(values) != 2:
            raise ValueError("OrderedPair expects format A~B")
        return OrderedPair(values[0], values[1])


def create_structure_db_info(structure):
    # Figure out the symmetry group
    sg = SymmetryFinder(structure, normalised_symmetry_precision(structure), -1)
    return create_structure_db_info_sg(structure, sg)


def create_structure_db_info_sg(structure, sg):
    d = dict()
    # Set the composition and formulas for the system
    comp = structure.composition
    el_amt = structure.composition.get_el_amt_dict()
    d.update({"unit_cell_formula": comp.to_dict,
              "reduced_cell_formula": comp.to_reduced_dict,
              "elements": list(el_amt.keys()),
              "nelements": len(el_amt),
              "pretty_formula": comp.reduced_formula,
              "anonymous_formula": comp.anonymized_formula,
              "nsites": comp.num_atoms,
              "chemsys": "-".join(sorted(el_amt.keys()))})

    d["spacegroup"] = {"symbol": unicode(sg.get_spacegroup_symbol(),
                                         errors="ignore"),
                       "number": sg.get_spacegroup_number(),
                       "point_group": unicode(sg.get_point_group(),
                                              errors="ignore"),
                       "source": "spglib",
                       "crystal_system": sg.get_crystal_system(),
                       "hall": sg.get_hall()}

    return d


def init_logging(args, logger):
    """Initialize verbosity
    """
    import logging

    logger.propagate = False
    hndlr = logging.StreamHandler()
    hndlr.setFormatter(logging.Formatter("[%(levelname)-6s] %(asctime)s %(name)s :: %(message)s"))
    logger.addHandler(hndlr)
    if args.quiet:
        lvl = logging.CRITICAL
    else:
        # Level:  default      -v            -vv
        lvl = (logging.WARN, logging.INFO, logging.DEBUG)[min(args.vb, 2)]
    logger.setLevel(lvl)

    # Log things from the splpy library
    splpy_logger = logging.getLogger("splpy")
    splpy_logger.addHandler(hndlr)
    splpy_logger.setLevel(lvl)


def normalised_symmetry_precision(structure, precision=0.01):
    len_per_site = (structure.volume / structure.num_sites) ** (1.0 / 3.0)
    dist_mtx = structure.distance_matrix
    np.fill_diagonal(dist_mtx, sys.float_info.max)
    shortest_dist = dist_mtx.min()
    if shortest_dist < len_per_site:
        return 0.1 * precision * shortest_dist

    return precision * len_per_site


def is_structure_bad(structure):
    # Check for max ratio between lattice parameters
    abc = structure.lattice.abc
    for i in range(0, 3):
        for j in range(i + 1, 3):
            if (max(abc[i], abc[j]) / min(abc[i], abc[j])) > 100:
                return True

    # Has the structure exploded?
    vol = structure.lattice.volume / len(structure)
    if vol > 100000:
        return True

    # Has the structure collapsed?
    lattice_mtx = structure.lattice.matrix
    for row in lattice_mtx:
        row /= LA.norm(row)
    norm_vol = np.dot(lattice_mtx[0], np.cross(lattice_mtx[1], lattice_mtx[2])) / len(structure)
    if norm_vol < 1e-2:
        return True

    return False


def write_structures(structures, names):
    str_list = [structures] if isinstance(structures, mg.Structure) else structures
    names_list = [names] if isinstance(names, str) else names

    for structure, name in zip(str_list, names_list):
        res = Res(structure)
        res.write_file("{}.res".format(name))
