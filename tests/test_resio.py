#!/usr/bin/env python

"""
Created on Apr 17, 2012
"""

from __future__ import division

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014, Martin Uhrin"
__version__ = "0.1"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "Feb 14, 2014"

import unittest
import os

from pymatgen import Composition, Lattice, Structure
from splpy.io.resio import Res

test_dir = os.path.dirname(__file__)

class ResTest(unittest.TestCase):

    def test_init(self):
        filepath = os.path.join(test_dir, 'input/23221-ZDsSsJoEW14.res')
        res = Res.from_file(filepath)
        comp = res.structure.composition
        self.assertEqual(comp, Composition.from_formula("C194H60"))
        #print res

        res_string = """TITL
CELL 1.0 1.0 1.0 1.0 90.0 90.0 90.0
LATT -1
SFAC Si F
Si 1 0.000000 0.000000 0.000000 1.0
F 2 0.750000 0.500000 0.750000 1.0"""
        res = Res.from_string(res_string)
        self.assertEqual(res.structure.composition, Composition("SiF"))
        self.assertEquals(res.structure.num_sites, 2)
        #print res

        struct = Structure(Lattice.orthorhombic(2.5, 3.5, 7.0), ['Na', 'Cl'], [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        res = Res(struct)
        res_string = str(res)
        lines = res_string.splitlines()
        self.assertEqual(lines[1], "CELL 1.0 2.5 3.5 7.0 90.0 90.0 90.0")


if __name__ == "__main__":
    unittest.main()
