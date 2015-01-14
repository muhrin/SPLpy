"""
File containing extensions to the pymatgen structure matcher
"""

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014, Martin Uhrin"
__version__ = "0.0.1"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "Jan 10, 2015"

import pymatgen.analysis.structure_matcher as structure_matcher
from pymatgen.symmetry.finder import SymmetryFinder

import splpy.util


class SymmetryComparator(structure_matcher.AbstractComparator):
    def __init__(self, parent=structure_matcher.SpeciesComparator()):
        self._parent = parent

    def are_equal(self, sp1, sp2):
        return self._parent.are_equal(sp1, sp2)

    def get_structure_hash(self, structure):
        """
        Hash for structure.

        Args:
            structure: A structure

        Returns:
            Reduced formula for the structure is used as a hash for the
            SpeciesComparator.
        """
        if not hasattr(structure, 'spacegroup'):
            sg = SymmetryFinder(structure, splpy.util.normalised_symmetry_precision(structure))
            structure.spacegroup = {"symbol": unicode(sg.get_spacegroup_symbol(), errors="ignore"),
                                    "number": sg.get_spacegroup_number(),
                                    "point_group": unicode(sg.get_point_group(), errors="ignore"),
                                    "crystal_system": sg.get_crystal_system(),
                                    "hall": sg.get_hall()}

        return "{} {}".format(structure.spacegroup['number'], self._parent.get_structure_hash(structure))