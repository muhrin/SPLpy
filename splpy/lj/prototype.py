"""
This module contains code relating to structure prototypes and their storage
in the database
"""
from splpy.util import normalised_symmetry_precision

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014"
__version__ = "0.1.0"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "Dec 22, 2014"

import itertools

from pymatgen.core import Structure
from pymatgen.transformations.standard_transformations import SubstitutionTransformation
from pymatgen.symmetry.finder import SymmetryFinder
import pymatgen.analysis.structure_matcher as structure_matcher

import splpy.util


def create_prototype(structure):
    reduced_comp = structure.composition.reduced_composition
    # Sort the composition from lowest number of particular species to highest
    els = sorted(reduced_comp.elements, key=lambda e: reduced_comp[e])

    # Convert the structure species to be A, B, C etc ordered in increasing quantity
    trans = SubstitutionTransformation(dict(zip(els, map(chr, range(ord('A'), ord('A') + len(els))))))
    str_copy = trans.apply_transformation(structure)
    # Scale the volume to be 1 unit per atom
    str_copy.scale_lattice(str_copy.num_sites)

    sg = SymmetryFinder(str_copy,
                        normalised_symmetry_precision(str_copy), -1)
    str_copy = sg.get_conventional_standard_structure()

    return str_copy


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
    prototypes = db['prototypes']

    structure = create_prototype(structure)
    # Find spacegroup
    sg = SymmetryFinder(structure,
                        normalised_symmetry_precision(structure), -1)

    matcher = structure_matcher.StructureMatcher(primitive_cell=False, attempt_supercell=True)
    transformations = create_transformations(structure)

    # Loop over all prototypes with correct stoichiometry and spacegroup
    for entry in prototypes.find(
            {"anonymous_formula": structure.composition.anonymized_formula,
             "spacegroup.number": sg.get_spacegroup_number()}):

        # Create the structure object for the prototype to compare against
        prototype = Structure.from_dict(entry["structure"])

        for transformation in transformations:
            structure_copy = transformation.apply_transformation(structure)

            # Compare with structure matcher (volume and species agnostic)
            if matcher.fit(prototype, structure_copy):
                return entry['_id']

    return None


def insert_prototype(structure, db):
    proto_id = find_prototype(structure, db)
    if proto_id:
        return proto_id, False

    prototype = create_prototype(structure)
    d = {"structure": prototype.to_dict}
    d.update(splpy.util.create_structure_db_info(prototype))
    proto_id = db['prototypes'].insert(d)

    return proto_id, True
