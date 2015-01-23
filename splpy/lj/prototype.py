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

COLLECTION_NAME = 'prototypes'

logger = logging.getLogger(__name__)


def create_prototype(structure):
    reduced_comp = structure.composition.reduced_composition
    # Sort the composition from lowest number of particular species to highest
    els = sorted(reduced_comp.elements, key=lambda e: reduced_comp[e])

    # Convert the structure species to be A, B, C etc ordered in increasing quantity
    trans = SubstitutionTransformation(dict(zip(els, map(chr, range(ord('A'), ord('A') + len(els))))))
    str_copy = trans.apply_transformation(structure)
    # Scale the volume to be 1 unit per atom
    str_copy.scale_lattice(str_copy.num_sites)

    sg = SymmetryFinder(str_copy, splpy.util.normalised_symmetry_precision(str_copy))

    return sg.get_conventional_standard_structure()


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
    # Find spacegroup
    sg = SymmetryFinder(structure, splpy.util.normalised_symmetry_precision(structure))

    matcher = structure_matcher.StructureMatcher(primitive_cell=False, attempt_supercell=True)
    transformed = [t.apply_transformation(structure) for t in create_transformations(structure)]

    # Loop over all prototypes with correct stoichiometry and spacegroup
    for entry in prototypes.find(
            {"anonymous_formula": structure.composition.anonymized_formula,
             "spacegroup.number": sg.get_spacegroup_number(),
             "nsites": structure.num_sites}):

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
    sg = SymmetryFinder(structure)
    d = {"structure": prototype.to_dict, "wyckoff_sites": get_wyckoff_sites(sg)}
    d.update(splpy.util.create_structure_db_info_sg(prototype, sg))
    proto_id = db[COLLECTION_NAME].insert(d)

    return proto_id, True

