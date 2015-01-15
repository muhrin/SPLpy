#!/usr/bin/env python

"""
Module implementing an Res file object class.
"""

from __future__ import division

__author__ = "Martin Uhrin, Georg Schusteritsch"
__copyright__ = "Copyright 2014, Martin Uhrin"
__version__ = "0.1.2"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "Feb 14, 2014"

import re

from monty.io import zopen

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.serializers.json_coders import MSONable


class Res(MSONable):
    """
    Object for representing the data in a Res file.
    Please note that this current implementation. Most attributes can be set
    directly.

    Args:
        structure (Structure):  Structure object.

    .. attribute:: structure

        Associated Structure.

    .. attribute:: name

        The name of the structure.

    .. attribute:: pressure

        The external pressure.

    .. attribute:: energy

        The internal energy of the structure.

    .. attribute:: spacegroup

        The space group of the structure.

    .. attribute:: times_found

        The number of times the structure was found.
    """

    def __init__(self, structure, name=None, pressure=None, energy=None, spacegroup=None, times_found=None):
        self.structure_ = structure
        self.name = name
        self.pressure = pressure
        self.energy = energy
        self.spacegroup = spacegroup
        self.times_found = times_found

    @property
    def structure(self):
        """
        Returns structure associated with this Res.
        """
        return self.structure_

    @staticmethod
    def from_file(filename):
        """
        Reads a Res from a file.

        Args:
            filename (str): File name containing Res data.

        Returns:
            Res object.
        """
        with zopen(filename, "r") as f:
            return Res.from_string(f.read())

    @staticmethod
    def parse_title(line):
        info = dict()

        tokens = line.split()
        num_tokens = len(tokens)
        # 1 = Name
        if num_tokens <= 1:
            return info
        info['name'] = tokens[1]
        # 2 = Pressure
        if num_tokens <= 2:
            return info
        info["pressure"] = float(tokens[2])
        # 3 = Volume
        # 4 = Internal energy
        if num_tokens <= 4:
            return info
        info["energy"] = float(tokens[4])
        # 5 = Spin density, 6 - Abs spin density
        # 7 = Space group OR num atoms (new format ONLY)
        idx = 7
        if tokens[idx][0] != '(':
            idx += 1

        if num_tokens <= idx:
            return info
        info["spacegroup"] = tokens[idx][1:len(tokens[idx]) - 1]
        # idx + 1 = n, idx + 2 = -
        # idx + 3 = times found
        if num_tokens <= idx + 3:
            return info
        info["times_found"] = int(tokens[idx + 3])

        return info

    @staticmethod
    def from_string(data):
        """
        Reads a Res from a string.

        Args:
            data (str): String containing Res data.

        Returns:
            Res object.
        """
        abc = []
        ang = []
        sp = []
        coords = []
        info = dict()
        coord_patt = re.compile(
            "(\w+)\s+([0-9]+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)"
        )
        lines = data.splitlines()
        line_no = 0
        while line_no < len(lines):
            line = lines[line_no]
            tokens = line.split()
            if tokens:
                if tokens[0] == 'TITL':
                    info = Res.parse_title(line)
                elif tokens[0] == 'CELL' and len(tokens) == 8:
                    abc = map(float, tokens[2:5])
                    ang = map(float, tokens[5:8])
                elif tokens[0] == 'SFAC':
                    for atom_line in lines[line_no:]:
                        if line.strip() == 'END':
                            break
                        else:
                            match = coord_patt.search(atom_line)
                            if match:
                                sp.append(match.group(1))  # 1-indexed
                                coords.append(map(float, match.groups()[2:5]))  # 0-indexed
                        line_no += 1  # Make sure the global is updated
            line_no += 1

        return Res(Structure(Lattice.from_lengths_and_angles(abc, ang), sp, coords), info.get('name'),
                   info.get('pressure'), info.get('energy'), info.get('spacegroup'), info.get('times_found'))

    def get_string(self, significant_figures=6):
        """
        Returns a string to be written as a Res file.

        Args:
            significant_figures (int): No. of significant figures to
                output all quantities. Defaults to 6.

        Returns:
            String representation of Res.
        """

        # Title line
        lines = ['TITL ' + self.print_title()]

        # Cell
        abc_ang = self.structure_.lattice.lengths_and_angles
        cell = ' '.join(map(str, abc_ang[0])) + ' ' + ' '.join(map(str, abc_ang[1]))
        lines.append("CELL 1.0 " + cell)

        # Latt
        lines.append('LATT -1')

        # Atoms
        species_types = self.structure_.symbol_set
        lines.append("SFAC " + ' '.join(map(str, species_types)))

        fmtstr = "{{}} {{}} {{:.{0}f}} {{:.{0}f}} {{:.{0}f}} 1.0".format(significant_figures)
        for site in self.structure_:
            symbol = site.specie.symbol
            coords = site.frac_coords
            lines.append(
                fmtstr.format(symbol, species_types.index(symbol) + 1, coords[0], coords[1],
                              coords[2]))
        lines.append('END')

        return "\n".join(lines)

    def __str__(self):
        """
        String representation of Res file.
        """
        return self.get_string()

    def write_file(self, filename, **kwargs):
        """
        Writes Res to a file. The supported kwargs are the same as those for
        the Res.get_string method and are passed through directly.
        """
        with open(filename, "w") as f:
            f.write(self.get_string(**kwargs) + "\n")

    @property
    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "structure": self.structure_.to_dict,
                "name": self.name,
                "pressure": self.pressure,
                "energy": self.energy,
                "spacegroup": self.spacegroup,
                "times_found": self.times_found}

    @classmethod
    def from_dict(cls, d):
        return Res(Structure.from_dict(d["structure"]), d.get('name'),
                   d.get('pressure'), d.get('energy'), d.get('spacegroup'), d.get('times_found'))

    def print_title(self):
        tokens = [self.name, self.pressure, self.structure.volume, self.energy, 0.0, 0.0, self.structure.num_sites]
        if self.spacegroup:
            tokens.append('(' + self.spacegroup + ')')
        else:
            tokens.append('(P1)')
        if self.times_found:
            tokens.append('n - ' + str(self.times_found))
        else:
            tokens.append('n - 1')

        return ' '.join(map(str, tokens))
