"""
This module defines a set of functions that make querying the Lennard Jones
database easier and groups common tasks.
"""
from splpy.util import OrderedPair

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014"
__version__ = "0.0.2"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "Jul 16, 2014"

import copy

import pymongo

import pymatgen as mg

import splpy.interval as interval
import splpy.lj.db_query
import splpy.lj.util as util
from splpy.lj.util import Criteriable
import splpy.util


class InteractionRange(Criteriable):
    """
    Object for representing a range of parameters in the Lennard Jones
    parameter space.  Each of the parameters can either be:
    Interval            - match parameter values in this interval
    None                - this parameter can take any value
    (string convertible)- match exactly this value
    """

    def __init__(self, epsilon=None, sigma=None, m=12, n=6, cut=2.5):
        self.epsilon = epsilon
        self.sigma = sigma
        self.m = m
        self.n = n
        self.cut = cut

    def to_criteria(self):
        criteria = dict()

        # Epsilon
        condition = self._param_to_condition(self.epsilon)
        if condition is not None:
            criteria["epsilon"] = condition

        # Sigma
        condition = self._param_to_condition(self.sigma)
        if condition is not None:
            criteria["sigma"] = condition

        # M
        condition = self._param_to_condition(self.m)
        if condition is not None:
            criteria["m"] = condition

        # N
        condition = self._param_to_condition(self.n)
        if condition is not None:
            criteria["n"] = condition

        # Cut
        condition = self._param_to_condition(self.cut)
        if condition is not None:
            criteria["cut"] = condition

        return criteria

    def _param_to_condition(self, param):
        if param is None:
            return None
        elif isinstance(param, interval.Interval):
            d = dict()
            if param.left_closed:
                d["$gte"] = param.start
            else:
                d["$gt"] = param.start

            if param.left_closed:
                d["$lte"] = param.end
            else:
                d["$lt"] = param.end

            return d
        else:
            return param

    @classmethod
    def from_dict(cls, d):
        parsed = dict()
        for val in ["epsilon", "sigma", "m", "n", "cut"]:
            if val in d:
                try:
                    parsed[val] = interval.Interval.from_string(d[val])
                except (ValueError, TypeError):
                    parsed[val] = float(d[val])

        return InteractionRange(**parsed)

    def __str__(self):
        conditions = list()
        if self.epsilon is not None:
            conditions.append("epsilon: {}".format(self.epsilon))
        if self.sigma is not None:
            conditions.append("sigma: {}".format(self.sigma))
        if self.m is not None:
            conditions.append("m: {}".format(self.m))
        if self.n is not None:
            conditions.append("n: {}".format(self.n))
        if self.cut is not None:
            conditions.append("cut: {}".format(self.cut))
        return ' '.join(conditions)


class LennardJonesSearchRange(object):

    def __init__(self):
        self.interactions = dict()

    def __str__(self):
        inters = ["{}: {{{}}}".format(pair, inter) for pair, inter in self.interactions.iteritems()]
        return ', '.join(inters)

    def add_interaction(self, species1, species2, inter):
        self.interactions[OrderedPair(species1, species2)] = inter

    def to_criteria(self):
        c = dict()
        for species_pair, crit in self.interactions.items():
            criteria = self._prepend_keys(crit.to_criteria(), species_pair)
            if criteria is not None:
                c.update(criteria)

        return c

    def _prepend_keys(self, criteria, pre):
        crit = dict()
        for key, value in criteria.iteritems():
            crit["{}.{}".format(pre, key)] = value
        return crit

    @property
    def num_interactions(self):
        return len(self.interactions)

    @classmethod
    def from_dict(cls, d):
        r = LennardJonesSearchRange()
        for pair, inter in d.iteritems():
            species = OrderedPair.from_string(pair)
            r.add_interaction(species.first, species.second, InteractionRange.from_dict(inter))
        return r


class VisitationEngine(object):
    class RunInstance:
        def __init__(self, visitor, properties=None):
            self.visitor = visitor
            self.stack_idx = 0
            self.properties = properties

    def __init__(self, query_engine):
        self._query_engine = query_engine
        self._generators = []
        self._run = None

    def add(self, generator):
        self._generators.append(generator)

    def run_queries(self, visitor, properties=None, criteria=None):
        self._run = self.RunInstance(visitor, properties)
        # Set the initial criteria that applies to all
        crit = criteria if criteria is not None else dict()
        self.callback(crit)
        self._run = None

    def callback(self, criteria):
        if self._run.stack_idx < len(self._generators):
            self._run.stack_idx += 1
            self._generators[self._run.stack_idx - 1](self._query_engine, criteria, self.callback)
            self._run.stack_idx -= 1
        else:
            self._run.visitor(self._query_engine.query(self._run.properties, criteria))


def visit_distinct_formulae(query_engine, criteria, callback):
    cur = query_engine.query(properties=["pretty_formula"], criteria=criteria)
    formulas = cur.distinct("pretty_formula")
    cur.close()
    for formula in formulas:
        crit = dict(criteria)
        crit["pretty_formula"] = formula
        callback(crit)


def visit_distinct_spacegroups(query_engine, criteria, callback):
    cur = query_engine.query(properties=["spacegroup.number"], criteria=criteria)
    spacegroups = cur.distinct("spacegroup.number")
    cur.close()
    for sg in spacegroups:
        crit = dict(criteria)
        crit["spacegroup.number"] = sg
        callback(crit)


class VisitParamPoints(object):
    def __init__(self, params=None):
        self._params_criteria = None
        if params:
            self._params_criteria = params.to_criteria()

    def __call__(self, query_engine, criteria, callback):
        for params_id in query_engine.get_param_ids(self._params_criteria):
            crit = copy.deepcopy(criteria) if criteria else dict()
            util.add_to_criteria(crit, "potential.params_id", params_id)
            callback(crit)

    def num_points(self, query_engine):
        return len(query_engine.get_param_ids(self._params_criteria))


def get_unique(query_engine, params, matcher, criteria=None, limit=None, save_doc=True):
    """
    Get the unique structures at a parameter point or points.  The params parameter can
    be either an LjInteractions or a InteractionRange.
    """
    class LowestEnergyStore:
        def __init__(self, lowest, matcher, limit):
            self.lowest = lowest
            self.matcher = matcher
            self.limit = limit

        def __call__(self, results):
            results.sort('energy', pymongo.ASCENDING)

            num_kept = 0
            for doc in results:
                formula = doc["pretty_formula"]
                spacegroup = doc["spacegroup.number"]

                formula_dict = self.lowest.setdefault(formula, dict())
                entries = formula_dict.setdefault(spacegroup, list())

                structure = mg.Structure.from_dict(doc["structure"])
                if not splpy.util.is_structure_bad(structure):
                    # Store the document with the structure so we can use it later
                    structure.splpy_doc = doc

                    # Check if we've seen this structure before
                    keep = True
                    for entry in entries:
                        try:
                            if matcher.fit(entry, structure):
                                keep = False
                                break
                        except MemoryError:
                            print("Ran out of memory trying to compare structures, can't get unique.")
                            print("Bad structures saved to local folder.")
                            splpy.util.write_structures(structure, [structure.splpy_doc["_id"] for structure in strs])
                            keep = False
                            break

                    if keep:
                        entries.append(structure)
                        num_kept += 1

                if self.limit and num_kept > self.limit:
                    break

    properties = ["_id", "name", "times_found", "energy", "spacegroup.symbol", "spacegroup.number", "pressure",
                  "structure", "pretty_formula"]

    crit = criteria if criteria else dict()
    ve = VisitationEngine(query_engine)
    ve.add(splpy.lj.db_query.VisitParamPoints(params))
    ve.add(splpy.lj.db_query.visit_distinct_formulae)

    lowest = dict()
    getter = LowestEnergyStore(lowest, matcher, limit)
    ve.run_queries(getter, properties=properties, criteria=crit)

    unique = list()
    while len(lowest) > 0:
        formula, spacegroups = lowest.popitem()
        for sg, structures in spacegroups.iteritems():
            unique.extend(structures)

    return unique


def surrounding_range(params, dist):
    """
    Get the InteractionRange that surround this point in parameter space up to a distance dist
    in each parameter direction
    """
    range = LennardJonesSearchRange()
    for pair, val in params.interactions.iteritems():
        d = val.__dict__
        intervals = dict()
        for param in ["epsilon", "sigma", "m", "n", "cut"]:
            val = d[param]
            intervals[param] = interval.Closed(val - dist, val + dist)
        range.add_interaction(pair.first, pair.second, InteractionRange(**intervals))

    return range