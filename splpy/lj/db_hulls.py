"""
This module contains code relating to convex hulls and their storage
in the database
"""

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014"
__version__ = "0.1.0"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "Jan 20, 2015"

import os
import StringIO
import shutil
import subprocess
import tempfile

import bson.objectid

import splpy.util as util
import splpy.lj as lj
import splpy.lj.db_structures as db_structures

HULLS_COLLECTION = 'hulls'


def get_stable_stoichiometries(db, params_id):
    ensure_hull(db, params_id)
    hulls_coll = db[HULLS_COLLECTION]
    return hulls_coll.find({"params_id": params_id, "entries.dist_above_hull": 0},
                           fields={"entries.pretty_formula"}).distinct("entries.pretty_formula")


def get_stable_structure_ids(db, params_id):
    ensure_hull(db, params_id)
    hulls_coll = db[HULLS_COLLECTION]
    return hulls_coll.find({"params_id": params_id, "entries.dist_above_hull": 0},
                           fields={"entries.structure_id"}).distinct("entries.structure_id")


def get_endpoints(db, params_id):
    structures_coll = db[db_structures.STRUCTURES_COLLECTION]
    # Build up the crtieria
    crit = db_structures.get_params_id_criteria(params_id)
    crit['elements'] = {'$size': 1}

    return structures_coll.find(crit).distinct("elements")


def ensure_hull(db, params_id):
    hulls = db[HULLS_COLLECTION]

    if not get_endpoints(db, params_id):
        hulls.remove({'params_id': params_id})
        return

    # Find or create the hull entry and find the id
    hull = util.find_or_create(hulls, {'params_id': params_id}, {'params_id': params_id},
                               fields={'last_updated': 1, 'entries.structure_id': 1})

    # Check if we need to remove any structures from the hull because they're no longer
    # in the structures collection
    structures = db['structures']
    structure_ids = [entry['structure_id'] for entry in hull['entries']] if 'entries' in hull else list()

    if len(structure_ids) != structures.find({'potential.params_id': params_id}).count() or \
       structures.find({'potential.params_id': params_id, '_id': {'$nin': structure_ids}}).limit(1).count() > 0:
        _regenerate_hull(db, hull['_id'])


def _regenerate_hull(db, hull_id):
    # Get all the structures at the parameter point of this hull
    hulls_coll = db[HULLS_COLLECTION]
    hull_doc = hulls_coll.find({'_id': hull_id})[0]
    params_id = hull_doc['params_id']

    endpoints = get_endpoints(db, params_id)
    if not endpoints:
        return None

    structures_coll = db[db_structures.STRUCTURES_COLLECTION]
    docs = list(structures_coll.find(db_structures.get_params_id_criteria(params_id)))

    # Save them to a temporary directory
    hull_dir = tempfile.mkdtemp()
    paths = list()
    for doc in docs:
        res = lj.util.create_writer(doc, "res")
        path = "{}.res".format(doc["_id"])
        res.write_file(os.path.join(hull_dir, path))
        paths.append(path)

    # Calculate formation enthalpies and distance above hull
    proc = subprocess.Popen(["sinfo", "-nf", "-i", "$f$,$rfo$,$hf$,$hd$\n"],
                            cwd=hull_dir, stdout=subprocess.PIPE, stdin=subprocess.PIPE, universal_newlines=True)
    output = StringIO.StringIO()
    output.write(proc.communicate(" ".join(paths))[0])
    # Delete the temporary folder
    shutil.rmtree(hull_dir)
    output.seek(0, 0)

    hull_entries = list()
    for line in output:
        tokens = line[:-1].split(',')
        entry = {"structure_id": bson.ObjectId(tokens[0][:-4]),
                 "pretty_formula": tokens[1],
                 "formation_energy": float(tokens[2]),
                 "dist_above_hull": float(tokens[3])}
        hull_entries.append(entry)

    # Save the new entries
    hulls_coll.find_and_modify({'_id': hull_id}, {"$set": {"endpoints": endpoints, "entries": hull_entries}},
                               upsert=True)