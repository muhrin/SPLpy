"""
This module defines a set of functions and classes that aid in manipulating
the Lennard-Jones structure database
"""
from splpy.util import normalised_symmetry_precision

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014"
__version__ = "0.0.1"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "Sept 23, 2014"

import copy
import datetime
import logging

from bson.objectid import ObjectId

from pymatgen.symmetry.finder import SymmetryFinder

import splpy.util

logger = logging.getLogger(__name__)


class LjDb:
    STRUCTURES_COLL = "structures"
    DUPLICATES_COLL = "duplicates"
    PARAMS_COLL = "params"
    PROTOTYPES_COLL = "prototypes"
    DUPLICATES_KEY = "file_name"


class ParamsInfo(object):
    def __init__(self, db, params):
        self._collection = db[LjDb.PARAMS_COLL]
        self.params = copy.copy(params)
        # Get rid of the MSON stuff
        self.params.pop("@module", None)
        self.params.pop("@class", None)
        self._params_id = params["_id"] if "_id" in params else None

    def fetch_id(self):
        if self._params_id is None:
            self._params_id = self._find_or_create_params_id()
        return self._params_id

    def _find_or_create_params_id(self):
        """
        Get the ObjectId for the parameters supplied.  If not found, create.
        """
        entry = splpy.util.find_or_create(self._collection, self.params, self.params, fields={"_id": 1})
        return entry["_id"]


def insert_structure(db, structure, params, name, energy, pressure):
    insert_structure_entry(db, _generate_entry(structure, params, name, energy, pressure), False)


def insert_structure_entry(db, entry, check_duplicates=True, update_duplicates=False):
    if "potential" not in entry or not "params" in entry["potential"]:
        return False

    entry["last_updated"] = datetime.datetime.today()

    # The collection we'll be inserting into
    coll = db[LjDb.STRUCTURES_COLL]

    params_info = ParamsInfo(db, entry["potential"]["params"])

    if check_duplicates:
        if LjDb.DUPLICATES_KEY not in entry:
            logger.error("Asked to check for duplicates on structure insert but duplicates key, {}, is missing".format(
                LjDb.DUPLICATES_KEY))
            return

        # WARNING: The version of the insert below is NOT atomic
        duplicates_value = entry[LjDb.DUPLICATES_KEY]
        result = coll.find_one({LjDb.DUPLICATES_KEY: duplicates_value}, fields=["_id"])

        if result is not None and not update_duplicates:
            logger.info("Skipping duplicate {} with id {}".format(duplicates_value, result["_id"]))
        else:
            # Need to either update or insert
            entry["potential"]["params_id"] = params_info.fetch_id()
            if result is None:
                entry["_id"] = ObjectId()
                coll.insert(entry)
            else:
                coll.update({"_id": result["_id"]}, {"$set": entry})
                # Make sure this comes after update or it will fail!
                entry.update(result)

            logger.info("Inserted {} with _id = {}".format(entry[LjDb.DUPLICATES_KEY], entry["_id"]))
            return entry["_id"]
    else:
        # Just insert - this is much faster
        entry["_id"] = ObjectId()
        entry["potential"]["params_id"] = params_info.fetch_id()
        coll.insert(entry)

        logger.debug("Inserted structure with id {}".format(entry["_id"]))
        return entry["_id"]


def _generate_entry(structure, params, name, energy, pressure):
    comp = structure.composition
    entry = {"structure": structure.to_dict, "name": name, "energy": energy, "energy_per_site": energy / comp.num_atoms,
             "pressure": pressure, "potential": {"name": "lennard_jones", "params": params.to_dict}}

    # Set the composition and formulas for the system
    el_amt = comp.get_el_amt_dict()
    entry.update({"unit_cell_formula": comp.to_dict,
                  "reduced_cell_formula": comp.to_reduced_dict,
                  "elements": list(el_amt.keys()),
                  "nelements": len(el_amt),
                  "pretty_formula": comp.reduced_formula,
                  "anonymous_formula": comp.anonymized_formula,
                  "nsites": comp.num_atoms,
                  "chemsys": "-".join(sorted(el_amt.keys()))})

    # Figure out the symmetry group
    sg = SymmetryFinder(structure, normalised_symmetry_precision(structure), -1)
    entry["spacegroup"] = {"symbol": str(sg.get_spacegroup_symbol(), errors="ignore"),
                           "number": sg.get_spacegroup_number(),
                           "point_group": str(sg.get_point_group(), errors="ignore"),
                           "source": "spglib",
                           "crystal_system": sg.get_crystal_system(),
                           "hall": sg.get_hall()}

    return entry

