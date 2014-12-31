#!/usr/bin/env python

"""
Created on Jul 22, 2014
"""

from __future__ import division
from splpy.util import OrderedPair

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014, Martin Uhrin"
__version__ = "0.1"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "Jul 22, 2014"

import unittest

from splpy.lj.db_query import LennardJonesSearchRange, HalfOpenInterval


class OrderedPairTest(unittest.TestCase):

    def test_init(self):
        ab = OrderedPair("A", "B")
        self.assertEqual(str(ab), "A~B")

class LennardJonesSearchRangeTest(unittest.TestCase):

    def test_init(self):

        r = LennardJonesSearchRange()

        interactions = list()
        interactions.append(OrderedPair("A", "B"))
        interactions.append(1.0)
        interactions.append(HalfOpenInterval(1.0, 1.5))
        interactions.append(12)
        interactions.append(6)
        interactions.append(2.5)

        r.add_interaction(LennardJonesSearchRange.InteractionRange(*interactions))

        print r.to_criteria()



if __name__ == "__main__":
    unittest.main()
