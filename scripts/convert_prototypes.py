

import argparse
import json
import sys

import bson.objectid

from pymatgen.core.structure import Structure
from pymatgen.symmetry.finder import SymmetryFinder

import matgendb as mgdb
import matgendb.util

from splpy.lj.query_engine import LjQueryEngine
import splpy.lj.prototype as prototype
import splpy.util


def parse_criteria_or_die(crit):
    criteria = None
    if crit:
        try:
            criteria = json.loads(crit)
        except ValueError:
            print(("Criteria {} is not a valid JSON string!".format(crit)))
            sys.exit(-1)
    return criteria

def init_criteria(args):
    criteria = parse_criteria_or_die(args.criteria)
    if criteria is None:
        criteria = dict()

    try:
        if args.stoichiometry:
            criteria['pretty_formula'] = args.stoichiometry.reduced_formula
    except AttributeError:
        pass
    try:
        if args.id:
            criteria['_id'] = bson.objectid.ObjectId(args.id)
    except AttributeError:
        pass

    for key, value in criteria.items():
        try:
            if value.startswith("ObjectId"):
                criteria[key] = bson.objectid.ObjectId(value[10:len(value) - 2])
        except AttributeError:
            pass

    return criteria

parser = argparse.ArgumentParser(description="""
Check prototypes wyckoff site match
""")

parser.add_argument('--quiet', '-q', dest='quiet', action="store_true", default=False,
                    help="Minimal verbosity.")
parser.add_argument('--verbose', '-v', dest='vb', action="count", default=0,
                    help="Print more verbose messages to standard error. Repeatable. (default=ERROR)")
parser.add_argument("-c", "--config", dest="config_file", type=str, default="ljdb.json",
                    help="Config file to use. Generate one using mgdb "
                         "init --config filename.json if necessary. "
                         "Otherwise, the code searches for an ljdb.json first"
                         "then ~/.ljdb.json. If none is found, an no-authentication "
                         "localhost:27017/lj database and structures "
                         "collection is assumed.")
parser.add_argument("--crit", dest="criteria", type=str, default=None, required=False,
                    help="Query criteria in typical json format.")

# Parse args
args = parser.parse_args()

d = mgdb.util.get_settings(args.config_file)
qe = LjQueryEngine(host=d["host"], port=d["port"], database=d["database"],
                   user=d["readonly_user"], password=d["readonly_password"],
                   collection=d["collection"], aliases_config=d.get("aliases_config", None))

db = qe.db
prototypes = db['prototypes']

criteria = init_criteria(args)
updated = 0
for doc in prototypes.find(criteria):
    print(("Updating {}: ({})".format(doc["_id"], updated)))
    new_proto = prototype.create_prototype(Structure.from_dict(doc["structure"]))

    sg = SymmetryFinder(new_proto, 0.01, angle_tolerance=-1)
    d = {"structure": new_proto.to_dict, "wyckoff_sites": prototype.get_wyckoff_sites(sg)}
    d.update(splpy.util.create_structure_db_info(new_proto, sg))
    prototypes.update({"_id": doc["_id"]}, {"$set": d})
    # splpy.util.write_structures([new_proto], ["proto"])

    updated += 1


