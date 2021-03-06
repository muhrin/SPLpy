"""
This module contains code relating to structures collection in the database
"""

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014"
__version__ = "0.1.0"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "Jan 23, 2015"

STRUCTURES_COLLECTION = 'structures'
DUPLICATES_COLLECTION = 'duplicates'


def get_params_id_criteria(params_id):
    return {'potential.params_id': params_id}


def get_params_id(structure_id, db):
    structures_coll = db[STRUCTURES_COLLECTION]
    cur = structures_coll.find({"_id": structure_id, "potential.params_id": {"$exists": True}},
                               fields={"potential.params_id"}).limit(1)
    if cur.count() == 0:
        return None

    return cur[0]["potential"]["params_id"]


def remove_structure(structure_id, db):
    # The collection we'll be inserting into
    structures = db[STRUCTURES_COLLECTION]
    duplicates = db[DUPLICATES_COLLECTION]

    structures.remove({'_id': structure_id})
    duplicates.remove({'duplicate_of': structure_id}, multi=True)