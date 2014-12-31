#!/usr/bin/env python

"""
Module for helper functions and other things that don't fit neatly elsewhere.
"""

from pymatgen.symmetry.finder import SymmetryFinder


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

    # Figure out the symmetry group
    sg = SymmetryFinder(structure, normalised_symmetry_precision(structure), -1)
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


def normalised_symmetry_precision(structure, precision=0.01):
    len_per_site = (structure.volume / structure.num_sites) ** 0.5
    return precision * len_per_site