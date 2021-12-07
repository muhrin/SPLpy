from pymatgen.core import Structure
from splpy.lj.query_engine import LjQueryEngine
from pymatgen.io.ase import AseAtomsAdaptor
import json

qe = LjQueryEngine(collection="prototypes")

formula = "A2B3"
for res in qe.query(criteria={"pretty_formula": formula}):
    i = 0
    filename = "{}_{}.json".format(formula, i)
    with open(filename, 'w') as f:
        json.dumps(res['structure'])
    i += 1
