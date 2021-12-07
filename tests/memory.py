__author__ = 'martin'


import glob
import os

import pymatgen.analysis.structure_matcher as structure_matcher

import splpy.util
import splpy.resio as resio
from splpy.lj.query_engine import LjQueryEngine
import splpy.lj.db_query as db_query
import splpy.lj.util as util

qe = LjQueryEngine(host="localhost", port=2701, user="martin", password="martin")

matcher = structure_matcher.StructureMatcher()


interactions = util.LjInteractions()
interactions.add_interaction(splpy.util.OrderedPair("A", "A"), util.LjInteraction(1, 1, 12, 6, 2.5))
interactions.add_interaction(splpy.util.OrderedPair("A", "B"), util.LjInteraction(2.8, 1, 12, 6, 2.5))
interactions.add_interaction(splpy.util.OrderedPair("B", "B"), util.LjInteraction(2.5, 1, 12, 6, 2.5))
params_range = db_query.surrounding_range(interactions, 0.3)

db_query.get_unique(qe, params_range, matcher)

strs = [resio.Res.from_file(res).structure for res in glob.glob(os.path.join('/home/martin/temp/crawl_fail/', '*.res'))]
for structure in strs:
    if splpy.util.is_structure_bad(structure):
        print('Bad structure')

groups = matcher.group_structures(strs)