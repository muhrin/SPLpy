
import argparse

import matgendb as mgdb
import matgendb.util

from splpy.lj.query_engine import LjQueryEngine

parser = argparse.ArgumentParser(description="""
Find the diamond structures in the db
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


for doc in qe.collection.find({"energy_per_site": {"$exists": False}},
                              fields=["_id", "energy", "nsites"]):
    qe.collection.update({"_id": doc["_id"]}, {"$set": {"energy_per_site": doc["energy"] / doc["nsites"]}})

duplicates = qe.db['duplicates']
for doc in duplicates.find({"energy_per_site": {"$exists": False}},
                                    fields=["_id", "energy", "nsites"]):
    duplicates.update({"_id": doc["_id"]}, {"$set": {"energy_per_site": doc["energy"] / doc["nsites"]}})