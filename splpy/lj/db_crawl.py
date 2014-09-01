"""
This module contains functions used for/when crawling the database.
"""

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014"
__version__ = "0.0.1"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "Aug 30, 2014"

import datetime
import itertools
import logging
import os
import yaml

from pymatgen.core.structure import Structure
import pymatgen.analysis.structure_matcher as structure_matcher

import splpy.lj.db_query as db_query
import splpy.lj.util

logger = logging.getLogger(__name__)


def get_next_param_to_crawl(query_engine, params_range=None):
    """
    This function looks for the oldest (longest time since last crawled) point in the
    parameter range (or all parameters if None).  First it looks for any parameters
    that have never been crawled (and hence have no crawl.last_crawled entry) after which
    it looks for the oldest point
    """
    params = query_engine.params

    # First try and find any parameters that have never been crawled
    criteria = {"crawl": {"$exists": False}}
    if params_range:
        criteria.update(params_range.to_criteria())
    result = params.find_and_modify(query=criteria, update={"$set": {"crawl.last_crawled": datetime.datetime.today()}},
                                    new=True)

    if result is None:
        # Find the oldest document to be crawled
        criteria = dict()
        if params_range:
            criteria.update(params_range.to_criteria())
        result = params.find_and_modify(query=criteria,
                                        update={"$set": {"crawl.last_crawled": datetime.datetime.today()}},
                                        sort={"crawl.last_crawled": 1}, new=True)

    return result


class Prune(object):
    def __init__(self):
        self.matcher = structure_matcher.StructureMatcher()

    def __call__(self, params, query_engine):
        # TODO: Check if there are any new structures since we last pruned

        # Save the fact that we pruned
        query_engine.params.update({"_id": params["_id"]}, {"$set": {"crawl.last_pruned": datetime.datetime.today()}})

        # Get all the structures at this parameter point
        structures = query_engine.query(criteria={"potential.params_id": params["_id"]})

        sorted_structures = sorted(structures, key=self._prune_hash)

        structures_coll = query_engine.collection
        duplicates_coll = query_engine.db["duplicates"]

        # Now check for any duplicates
        for k, g in itertools.groupby(sorted_structures, key=self._prune_hash):
            unmatched = list()
            for s in g:
                structure = Structure.from_dict(s["structure"])
                structure.splpy_doc = s
                unmatched.append(structure)

            while len(unmatched) > 0:
                ref = unmatched.pop()

                # Collect together the matches
                inds = filter(lambda i: self._fit(ref, unmatched[i]), xrange(len(unmatched)))
                duplicates = [unmatched[i].splpy_doc for i in inds]
                # Save the id so we know what it's a duplicate of
                for duplicate in duplicates:
                    duplicate["duplicate_of"] = ref.splpy_doc["_id"]

                if len(duplicates) > 0:
                    logger.info(
                        "Found {} duplicates of structure with id {}".format(len(duplicates), ref.splpy_doc["_id"]))
                    structures_coll.remove({"_id": {"$in": [s["_id"] for s in duplicates]}})
                    duplicates_coll.insert(duplicates)

                    unmatched = [unmatched[i] for i in xrange(len(unmatched)) if i not in inds]

    def _prune_hash(self, structure_doc):
        return structure_doc["pretty_formula"], structure_doc["spacegroup"]["symbol"]

    def _fit(self, structure1, structure2, energy_tolerance=1e-2):

        doc1 = structure1.splpy_doc
        doc2 = structure2.splpy_doc
        if doc1["spacegroup"] != doc2["spacegroup"]:
            return False
        # Check the two energies are within the fractional energy tolerance
        if not are_close_fraction(doc1["energy"], doc2["energy"], energy_tolerance):
            return False

        return self.matcher.fit(structure1, structure2)


class Refine(object):
    def __init__(self, spipe_input, dist, limit=None):
        self.spipe_input = spipe_input
        self.dist = dist
        self.limit = limit
        # TODO: Check if spipe is installed

    def __call__(self, params_doc, query_engine):
        params = splpy.lj.util.LjInteractions.from_dict(params_doc)
        range = db_query.surrounding_range(params, self.dist)

        matcher = structure_matcher.StructureMatcher()
        strs = db_query.get_unique_in_range(query_engine, range, matcher, limit=self.limit)

        params_id = str(params_doc["_id"])
        try:
            run_dir = params_id
            structures_dir = os.path.join(run_dir, "unique")
            os.makedirs(structures_dir)

            for struct in strs:
                doc = struct.splpy_doc
                res = splpy.resio.Res(struct, doc.get("name"), doc.get("pressure"), doc.get("energy"),
                                      doc.get("spacegroup.symbol"), doc.get("times_found"))
                res.write_file(os.path.join(structures_dir, "{}.res".format(doc["_id"])))

            # Write the spipe yaml file for the run
            with open(os.path.join(run_dir, self.spipe_input), 'w') as outfile:
                spipe_settings = dict()
                if os.path.exists(self.spipe_input):
                    with open(self.spipe_input, 'r') as original:
                        spipe_settings = yaml.load(original)
                self._update_spipe_dict(params, "unique", spipe_settings)
                outfile.write(yaml.dump(spipe_settings))

                # TODO: Run spipe

                # TODO: Read in all generated structures

                # TODO: Find any new structures

                # TODO: Insert any new structures into the database


        except OSError:
            pass

    def _update_spipe_dict(self, params, structures_dir, d):
        d["loadStructures"] = structures_dir

        params_dict = dict()
        for pair, inter in params.interactions.iteritems():
            params_dict[str(pair)] = [inter.epsilon, inter.sigma, inter.m, inter.n, inter.cut]
        if not d.get("potential"):
            d["potential"] = dict()
        d["potential"]["lennardJones"] = {"params": params_dict}

        return d


def are_close_fraction(value1, value2, tolerance):
    return (value1 - value2) / value1 < tolerance and (value1 - value2) / value2 < tolerance