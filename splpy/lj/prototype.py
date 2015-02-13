"""
This module contains code relating to structure prototypes and their storage
in the database
"""

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014"
__version__ = "0.1.0"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "Dec 22, 2014"

import itertools
import logging

from pymatgen.core import Structure
from pymatgen.transformations.standard_transformations import SubstitutionTransformation
from pymatgen.symmetry.finder import SymmetryFinder
import pymatgen.analysis.structure_matcher as structure_matcher

import splpy.util
from splpy.util import normalised_symmetry_precision
import splpy.resio
import splpy.structure_tidy as structure_tidy

COLLECTION_NAME = 'prototypes'

logger = logging.getLogger(__name__)

_symm_precision = 0.01


def create_prototype(structure):
    reduced_comp = structure.composition.reduced_composition
    # Sort the composition from lowest number of particular species to highest
    els = sorted(reduced_comp.elements, key=lambda e: reduced_comp[e])

    # Convert the structure species to be A, B, C etc ordered in increasing quantity
    trans = SubstitutionTransformation(dict(zip(els, map(chr, range(ord('A'), ord('A') + len(els))))))
    str_copy = trans.apply_transformation(structure)
    # Scale the volume to be 1 unit per atom
    str_copy.scale_lattice(str_copy.num_sites)

    sg = SymmetryFinder(str_copy, _symm_precision, angle_tolerance=-1)
    # return structure_tidy.structure_tidy(sg.get_primitive_standard_structure())
    #return structure_tidy.structure_tidy(sg.find_primitive())
    #return structure_tidy.structure_tidy(sg.get_refined_structure())
    return sg.get_refined_structure()


def create_transformations(structure):
    reduced_comp = structure.composition.reduced_composition
    # Sort the composition from lowest number of particular species to highest
    els = sorted(reduced_comp.elements, key=lambda e: reduced_comp[e])

    # Split the composition into groups of the same number of species
    specie_groups = [[]]
    num_of_last_specie = reduced_comp[els[0]]
    for el in els:
        if num_of_last_specie == reduced_comp[el]:
            specie_groups[-1].append(el)
        else:
            specie_groups.append([el])

    mappings = [dict()]
    last_el = ord('A')
    for group in specie_groups:
        new_mappings = list()
        for perm in itertools.permutations(map(chr, range(last_el, last_el + len(group)))):
            for mapping in mappings:
                new_mappings.append(mapping.copy())
                new_mappings[-1].update(zip(group, perm))

        mappings = new_mappings
        last_el += len(group)

    return [SubstitutionTransformation(m) for m in mappings]


def find_prototype(structure, db):
    # Get the prototypes collection
    prototypes = db[COLLECTION_NAME]

    structure = create_prototype(structure)
    if not structure:
        return None
    # Find spacegroup
    sg = SymmetryFinder(structure, _symm_precision, angle_tolerance=-1)

    matcher = structure_matcher.StructureMatcher(primitive_cell=False, attempt_supercell=True)
    transformed = [t.apply_transformation(structure) for t in create_transformations(structure)]

    # Loop over all prototypes with correct stoichiometry and spacegroup
    for entry in prototypes.find(
            {"anonymous_formula": structure.composition.anonymized_formula,
             "spacegroup.number": sg.get_spacegroup_number(),
             "nsites": structure.num_sites,
             "wyckoff_sites": get_wyckoff_sites(sg)}):

        # Create the structure object for the prototype to compare against
        prototype = Structure.from_dict(entry["structure"])

        for trans in transformed:
            # Compare with structure matcher (volume and species agnostic)
            if matcher.fit(prototype, trans):
                return entry['_id']

    return None


def get_wyckoff_sites(spacegroup):
    class Counter(dict):
        def __missing__(self, key):
            return 0

    wyckoff_sites = Counter()
    for site in spacegroup.get_symmetry_dataset()['wyckoffs']:
        wyckoff_sites[site] += 1
    return wyckoff_sites


def insert_prototype(structure, db):
    proto_id = find_prototype(structure, db)
    if proto_id:
        return proto_id, False

    prototype = create_prototype(structure)
    sg = SymmetryFinder(prototype, _symm_precision, angle_tolerance=-1)
    d = {"structure": prototype.to_dict, "wyckoff_sites": get_wyckoff_sites(sg)}
    d.update(splpy.util.create_structure_db_info(prototype, sg))
    proto_id = db[COLLECTION_NAME].insert(d)

    return proto_id, True


def get_all_same(proto_id, db):
    # Get the prototypes collection
    prototypes = db[COLLECTION_NAME]

    cur = prototypes.find({'_id': proto_id,
                           '$or': [{'structure_type': {'$exists': True}}, {'strukturbericht': {'$exists': True}}]},
                          fields={'_id': False, 'structure_type': True, 'pearson': True, 'strukturbericht': True})
    if cur.count() == 0:
        return [proto_id]

    info = cur[0]
    return prototypes.find(info).distinct('_id')


def get_label(proto_id, db):
    # Get the prototypes collection
    prototypes = db[COLLECTION_NAME]

    cur = prototypes.find({'_id': proto_id},
                          fields={'_id': False, 'structure_type': True, 'pearson': True,
                                  'strukturbericht': True, 'extra_labels': True})
    if cur.count() == 0:
        return None

    d = cur[0]
    if 'strukturbericht' in d and d['strukturbericht']:
        return d['strukturbericht']
    elif 'structure_type' in d and d['structure_type']:
        return d['structure_type']
    elif 'extra_labels' in d and d['extra_labels']:
        return ", ".join(
            ["{}: {}".format(name, value) for name, value in d['extra_labels'].iteritems()])
    else:
        return str(proto_id)