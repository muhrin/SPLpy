"""
This module contains functions used for/when crawling the database.
"""

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014"
__version__ = "0.1.0"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "Aug 30, 2014"

import datetime
import glob
import itertools
import logging
import os
import subprocess
import shutil
import tempfile
import yaml

import pymongo

from pymatgen.core.structure import Structure
import pymatgen.analysis.structure_matcher as structure_matcher

import splpy.lj.db_hulls as db_hulls
import splpy.lj.db_manip as db_manip
import splpy.lj.db_query as db_query
import splpy.lj.db_structures as db_structures
import splpy.lj.prototype as prototype
import splpy.lj.util
import splpy.resio as resio
import splpy.util as util
import splpy.structure_matching
from splpy.lj.db_query import VisitationEngine

logger = logging.getLogger(__name__)


def get_next_param_to_crawl(query_engine, params_range=None, criteria=None):
    """
    This function looks for the oldest (longest time since last crawled) point in the
    parameter range (or all parameters if None).  First it looks for any parameters
    that have never been crawled (and hence have no crawl.last_crawled entry) after which
    it looks for the oldest point
    """
    params = query_engine.params

    crit = dict()
    if criteria:
        crit.update(criteria)

    # First try and find any parameters that have never been crawled
    crit.update({"crawl": {"$exists": False}})
    if params_range:
        crit.update(params_range.to_criteria())
    result = params.find_and_modify(query=crit, update={"$set": {"crawl.last_crawled": datetime.datetime.today()}},
                                    new=True)

    if result is None:
        # Find the oldest document to be crawled
        crit = dict()
        if criteria:
            crit.update(criteria)
        if params_range:
            crit.update(params_range.to_criteria())
        result = params.find_and_modify(query=crit,
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
            for doc in g:
                structure = Structure.from_dict(doc["structure"])
                structure.splpy_doc = doc
                if self._bad_structure(structure):
                    db_structures.remove_structure(doc["_id"], query_engine.db)
                else:
                    unmatched.append(structure)

            while len(unmatched) > 0:
                ref = unmatched.pop()

                # Collect together the matches
                inds = [i for i in range(len(unmatched)) if self._fit(ref, unmatched[i])]
                duplicates = [unmatched[i].splpy_doc for i in inds]
                # Save the id so we know what it's a duplicate of
                for duplicate in duplicates:
                    duplicate["duplicate_of"] = ref.splpy_doc["_id"]

                if duplicates:
                    logger.info(
                        "Found {} duplicates of structure with id {}".format(len(duplicates), ref.splpy_doc["_id"]))
                    structures_coll.remove({"_id": {"$in": [s["_id"] for s in duplicates]}})
                    duplicates_coll.insert(duplicates)

                    unmatched = [unmatched[i] for i in range(len(unmatched)) if i not in inds]

    def _prune_hash(self, structure_doc):
        return structure_doc["pretty_formula"], structure_doc["spacegroup"]["number"]

    def _bad_structure(self, structure):
        # Kill any structures that have a near zero energy
        if "energy" in structure.splpy_doc and structure.splpy_doc["energy"] > -0.01:
            return True
        if util.is_structure_bad(structure):
            return True
        return False

    def _fit(self, structure1, structure2, energy_tolerance=1e-2):

        doc1 = structure1.splpy_doc
        doc2 = structure2.splpy_doc
        if doc1["spacegroup"]["number"] != doc2["spacegroup"]["number"]:
            return False
        # Check the two energies are within the fractional energy tolerance
        if not are_close_fraction(doc1["energy"] / len(structure1), doc2["energy"] / len(structure2),
                                  energy_tolerance):
            return False

        return self.matcher.fit(structure1, structure2)


class Refine(object):
    def __init__(self, spipe_input, dist, limit=None):
        self.spipe_input = spipe_input
        self.dist = dist
        self.limit = limit

        if not os.path.exists(self.spipe_input):
            print(("Error: {} does not exist.".format(self.spipe_input)))

        # Check if spipe is installed
        try:
            from subprocess import DEVNULL  # py3k
        except ImportError:
            DEVNULL = open(os.devnull, 'wb')

        try:
            subprocess.call(["spipe", "--help"], stdout=DEVNULL)
            self.found_spipe = True
        except OSError:
            logger.error("spipe not found, can't refine database.")
            self.found_spipe = False

        # self._matcher = structure_matcher.StructureMatcher(comparator=splpy.structure_matching.SymmetryComparator())
        self._matcher = structure_matcher.StructureMatcher(primitive_cell=False)

    def __call__(self, params_doc, query_engine):
        if not self.found_spipe:
            return

        params_id = params_doc["_id"]
        params = splpy.lj.util.LjInteractions.from_dict(params_doc)

        # Get the uniques at the current parameter point
        my_structures = db_query.get_unique(query_engine, params, self._matcher,
                                            save_doc=True, properties=['prototype_id'])
        proto_ids = set()
        for structure in my_structures:
            proto_id = structure.splpy_doc.get("prototype_id")
            if proto_id:
                proto_ids.add(proto_id)

        # Get uniques from surrounding points (exclude this point)
        params_range = db_query.surrounding_range(params, self.dist)
        surrounding_structures = self._get_unique(params_range, query_engine, proto_ids,
                                                  criteria={"potential.params_id": {"$ne": params_id}})

        logger.debug("Found {} unique surrounding structures".format(len(surrounding_structures)))

        # Cull uniques that are the same as any of my structures
        for my in my_structures:
            for surrounding in surrounding_structures:
                if util.is_structure_bad(surrounding) or self._matcher.fit(my, surrounding):
                    surrounding_structures.remove(surrounding)
                    break

        logger.debug("{} unique surrounding structure(s) left after excluding those at this param point".
                     format(len(surrounding_structures)))

        # Run spipe on remaining unique structures
        try:
            run_dir = tempfile.mkdtemp()
        except OSError:
            logger.error("Failed to create temporary directory to run spipe")
            return

        spipe_input = self._prepare_spipe_files(params, surrounding_structures, run_dir)
        del surrounding_structures[:]

        # Run spipe
        subprocess.Popen(["spipe", spipe_input], cwd=run_dir).wait()

        # Cull any structures that have relaxed to one of mine

        # Read in all generated structures and group uniques
        # relaxed_structures = [resio.Res.from_file(res).structure for res in glob.glob(os.path.join(run_dir, '*.res'))]
        relaxed_structures = list()
        for file in glob.glob(os.path.join(run_dir, '*.res')):
            res = resio.Res.from_file(file)
            if not util.is_structure_bad(res.structure):
                res.structure.splpy_res = res
                relaxed_structures.append(res.structure)

        new_structures = list()
        while relaxed_structures:
            relaxed = relaxed_structures.pop()
            energy_per_atom = relaxed.splpy_res.energy / len(relaxed)

            keep = True
            # Check if it's the same as any of the structures at this parameter point already
            for my in my_structures:
                if are_close_fraction(energy_per_atom, my.splpy_doc['energy'] / len(my), 1e-2):
                    try:
                        if self._matcher.fit(relaxed, my):
                            keep = False
                            break
                    except MemoryError:
                        # Sometimes matcher runs out of memory if it tried to handle a structure that has
                        # many nearest neighbour interactions (e.g. a skewed cell)
                        logger.error("Memory error trying to match structures.")
                        keep = False

            if keep:
                # Check if it's the same as any of the other relaxed structures
                for unique in new_structures:
                    if are_close_fraction(energy_per_atom, unique.splpy_res.energy / len(unique), 1e-2):
                        try:
                            if self._matcher.fit(relaxed, unique):
                                keep = False
                                break
                        except MemoryError:
                            # Sometimes matcher runs out of memory if it tried to handle a structure that has
                            # many nearest neighbour interactions (e.g. a skewed cell)
                            logger.error("Memory error trying to match structures.")
                            keep = False

            if keep:
                new_structures.append(relaxed)

        logger.info("Found {} new structures at parameter point {}".format(len(new_structures), params_id))

        # Save remainder to db
        for new in new_structures:
            res = new.splpy_res
            db_manip.insert_structure(query_engine.db, new, params, res.name, res.energy, res.pressure)

        # Delete the temporary folder
        shutil.rmtree(run_dir)

    def _prepare_spipe_files(self, params, structures, dir):
        structures_dir = os.path.join(dir, "unique")
        os.makedirs(structures_dir)

        # Save the res files
        for struct in structures:
            doc = struct.splpy_doc
            res = resio.Res(struct, doc.get("_id"), doc.get("pressure"), doc.get("energy"),
                            doc.get("spacegroup.symbol"), doc.get("times_found"))
            res.write_file(os.path.join(structures_dir, "{}.res".format(doc["_id"])))

        # Write the spipe yaml file for the run
        spipe_input = os.path.basename(self.spipe_input)
        with open(os.path.join(dir, spipe_input), 'w') as outfile:
            spipe_settings = dict()
            if os.path.exists(self.spipe_input):
                with open(self.spipe_input, 'r') as original:
                    spipe_settings = yaml.load(original)
            self._update_spipe_dict(params, "unique", spipe_settings)
            outfile.write(yaml.dump(spipe_settings))

        return spipe_input

    def _update_spipe_dict(self, params, structures_dir, d):
        d["loadStructures"] = structures_dir

        params_dict = dict()
        for pair, inter in list(params.interactions.items()):
            params_dict[str(pair)] = [inter.epsilon, inter.sigma, inter.m, inter.n, inter.cut]
        potential = d.setdefault("potential", dict())
        potential["lennardJones"] = {"params": params_dict}

        return d

    def _get_unique(self, params, query_engine, known_prototypes, criteria=None):
        class LowestEnergyStore:
            def __init__(self, lowest, known_prototypes, limit):
                self.lowest = lowest
                self.limit = limit
                self.known_prototypes = known_prototypes

            def __call__(self, results):
                results.sort('energy_per_site', pymongo.ASCENDING)
                num_kept = 0
                for doc in results:
                    formula = doc["pretty_formula"]
                    formula_dict = self.lowest.setdefault(formula, list())

                    keep = False
                    proto_id = doc.get("prototype_id")
                    if proto_id:
                        if proto_id not in self.known_prototypes:
                            keep = True
                            self.known_prototypes.add(proto_id)
                    else:
                        keep = True

                    if keep:
                        structure = Structure.from_dict(doc["structure"])
                        structure.splpy_doc = doc
                        formula_dict.append(structure)
                        num_kept += 1

                    if self.limit and num_kept >= self.limit:
                        break

        crit = criteria if criteria else dict()

        ve = VisitationEngine(query_engine)
        params_query = splpy.lj.db_query.VisitParamPoints(params)
        logger.debug("Refining from {} surrounding points".format(params_query.num_points(query_engine)))
        ve.add(params_query)
        ve.add(splpy.lj.db_query.visit_distinct_formulae)

        structures = dict()
        getter = LowestEnergyStore(structures, known_prototypes, self.limit)
        ve.run_queries(getter, properties=["_id", "pretty_formula", "potential.params_id", "structure", "prototype_id"],
                       criteria=crit)

        all_structures = list()
        for formula, strs in list(structures.items()):
            all_structures.extend(strs)
        return all_structures


class AssignPrototypes(object):
    def __init__(self, limit):
        self.limit = limit

    def __call__(self, params, query_engine):
        # Assign prototypes to all the structures at this parameter point that do not have a prototype
        total_prototypes = 0
        new_prototypes = 0

        structures_coll = query_engine.collection
        stoichs = structures_coll.find({"potential.params_id": params["_id"]}, {"pretty_formula": 1}).\
            distinct("pretty_formula")

        for stoich in stoichs:
            cur = structures_coll.find({"potential.params_id": params["_id"], "prototype_id": {"$exists": False},
                                        "pretty_formula": stoich}, fields={"structure": 1, "energy_per_site": 1})
            if self.limit:
                cur.sort('energy_per_site', pymongo.ASCENDING).limit(self.limit)

            for doc in cur:
                # Either get the prototype or insert this structure as a new one
                structure = Structure.from_dict(doc["structure"])
                if not splpy.util.is_structure_bad(structure):
                    proto_id, is_new = prototype.insert_prototype(structure, query_engine.db)
                    if proto_id is None:
                        logger.warn("Failed to create find or create prototype for structure {}".format(doc["_id"]))
                        continue

                    query_engine.collection.update({'_id': doc['_id']}, {"$set": {'prototype_id': proto_id}})
                    total_prototypes += 1
                    if is_new:
                        new_prototypes += 1
                else:
                    logger.info("Skipping structure {} because there is something wrong with it.".format(doc["_id"]))

        logger.info("Assigned {} prototypes, {} new.".format(total_prototypes, new_prototypes))


def ensure_hull(params_doc, query_engine):
    db_hulls.ensure_hull(query_engine.db, params_doc['_id'])


def are_close_fraction(value1, value2, tolerance):
    # Guard against one or both value 1/2 being zero (because we divide later on)
    if value1 == 0 and value2 == 0:
        return True
    else:
        if value1 == 0 or value2 == 0:
            return False
        return abs(value1 - value2) / abs(value1) < tolerance and \
               abs(value1 - value2) / abs(value2) < tolerance
