import argparse

import pymongo

import pymatgen.symmetry.finder
import pymatgen.analysis.structure_matcher as structure_matcher
import pymatgen as mg

import matgendb as mgdb
import matgendb.util

from splpy.lj.query_engine import LjQueryEngine
import splpy.resio
from splpy.lj.db_query import VisitationEngine
import splpy.lj.util
from splpy.util import normalised_symmetry_precision


class DiamondFinder:
    def __init__(self):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = mg.Lattice.from_parameters(a=3.84, b=3.84, c=3.84, alpha=120,
                                             beta=90, gamma=60)
        self.diamond = mg.Structure(lattice, ["B", "B"], coords)

        self.matcher = structure_matcher.StructureMatcher()

    def __call__(self, results):
        # This is a bid dodgy (using private function), but QueryResults can't handle
        # calling through to Cursor methods properly (and not loose its identity)
        if len(results) == 0:
            return

        doc = results._mapped_result(results.sort('energy', pymongo.ASCENDING).limit(1)[0])

        structure_dict = doc['structure']
        for site in structure_dict['sites']:
            site['label'] = 'B'
            site['species'] = [{'occu': 1, 'element': 'B'}]
        structure = mg.Structure.from_dict(structure_dict)
        sg = mg.symmetry.finder.SymmetryFinder(structure, normalised_symmetry_precision(structure), -1)

        if sg.get_spacegroup_number() == 227 and self.matcher.fit(self.diamond, structure):
            params = doc['potential']['params']
            print(('Found diamond: {}'.format(doc['_id'])))
            print(('{} {} {} {}'.format(params['A~B']['epsilon'], params['A~B']['sigma'], params['B~B']['epsilon'],
                                       params['B~B']['sigma'])))
            res = splpy.resio.Res(structure, doc['_id'])
            res.write_file('{}.res'.format(doc['_id']))

        results.close()


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

ve = VisitationEngine(qe)
params_query = splpy.lj.db_query.VisitParamPoints()
ve.add(params_query)

finder = DiamondFinder()
ve.run_queries(finder, properties={"_id": 1, "energy": 1, "structure": 1, "potential": 1}, criteria=dict())
