

import argparse

import bson.objectid

from pymatgen.core.structure import Structure
from pymatgen.symmetry.finder import SymmetryFinder

import matgendb as mgdb
import matgendb.util

from splpy.lj.query_engine import LjQueryEngine
import splpy.lj.prototype as prototype
import splpy.structure_tidy as structure_tidy

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

# Parse args
args = parser.parse_args()

d = mgdb.util.get_settings(args.config_file)
qe = LjQueryEngine(host=d["host"], port=d["port"], database=d["database"],
                   user=d["readonly_user"], password=d["readonly_password"],
                   collection=d["collection"], aliases_config=d.get("aliases_config", None))

db = qe.db
structures = db['structures']
prototypes = db['prototypes']


str_doc = structures.find({"_id": bson.objectid.ObjectId("54e216cb2d4cd80d30afd6ed")})[0]
struct = Structure.from_dict(str_doc["structure"])
proto = prototype.create_prototype(struct)


precision = 0.01
same = 0
different = 0
for doc in structures.find({"prototype_id": {"$exists": True}, "spacegroup.number": {"$gt": 100}},
                           fields={"prototype_id": 1, "structure": 1}):
    proto_doc = prototypes.find({"_id": doc["prototype_id"]}, fields={"structure": 1, "wyckoff_sites": 1})[0]

    proto_pretidy = Structure.from_dict(proto_doc["structure"])
    sg = SymmetryFinder(proto_pretidy, precision, -1)
    proto = structure_tidy.structure_tidy(sg.get_primitive_standard_structure())
    sg_proto = SymmetryFinder(proto, precision, -1)
    proto_wyckoff = prototype.get_wyckoff_sites(sg_proto)

    str_pretidy = prototype.create_prototype(Structure.from_dict(doc["structure"]))
    sg = SymmetryFinder(str_pretidy, precision, -1)
    str = structure_tidy.structure_tidy(sg.get_primitive_standard_structure())
    sg_str = SymmetryFinder(str, precision, -1)
    str_wyckoff = prototype.get_wyckoff_sites(sg_str)
    # if str_wyckoff != proto_wyckoff:
    #     print("Wyckoff sites not equal!")
    #     print("proto: {}".format(proto_wyckoff))
    #     print("structure: {}".format(str_wyckoff))
    #     splpy.util.write_structures([proto, str], ['proto', 'str'])
    if str_wyckoff != proto_wyckoff:
        different += 1
        print("{} {}".format(str_wyckoff == proto_wyckoff, sg_str.get_spacegroup_number()))
    else:
        same += 1
    print("Same: {}, Different: {}, Percent: {}".format(same, different, float(different)/float(same+different)))
