#!/usr/bin/env python

"""
An interface for using structure tidy to standardise a structure

"""

from __future__ import division

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2012, Martin Uhrin"
__version__ = "0.1"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin@gmail.com"
__date__ = "Jan 27, 2015"

import os
import re
import subprocess
import StringIO
import sys
import tempfile

from pymatgen.core.structure import Structure
import pymatgen.core.lattice


# Copied from shutils Python 3.3
def which(cmd, mode=os.F_OK | os.X_OK, path=None):
    """Given a command, mode, and a PATH string, return the path which
    conforms to the given mode on the PATH, or None if there is no such
    file.

    `mode` defaults to os.F_OK | os.X_OK. `path` defaults to the result
    of os.environ.get("PATH"), or can be overridden with a custom search
    path.

    """
    # Check that a given file can be accessed with the correct mode.
    # Additionally check that `file` is not a directory, as on Windows
    # directories pass the os.access check.
    def _access_check(fn, mode):
        return (os.path.exists(fn) and os.access(fn, mode)
                and not os.path.isdir(fn))

    # If we're given a path with a directory part, look it up directly rather
    # than referring to PATH directories. This includes checking relative to the
    # current directory, e.g. ./script
    if os.path.dirname(cmd):
        if _access_check(cmd, mode):
            return cmd
        return None

    if path is None:
        path = os.environ.get("PATH", os.defpath)
    if not path:
        return None
    path = path.split(os.pathsep)

    if sys.platform == "win32":
        # The current directory takes precedence on Windows.
        if not os.curdir in path:
            path.insert(0, os.curdir)

        # PATHEXT is necessary to check on Windows.
        pathext = os.environ.get("PATHEXT", "").split(os.pathsep)
        # See if the given file matches any of the expected path extensions.
        # This will allow us to short circuit when given "python.exe".
        # If it does match, only test that one, otherwise we have to try
        # others.
        if any(cmd.lower().endswith(ext.lower()) for ext in pathext):
            files = [cmd]
        else:
            files = [cmd + ext for ext in pathext]
    else:
        # On other platforms you don't have things like PATHEXT to tell you
        # what file suffixes are executable, so just pass on cmd as-is.
        files = [cmd]

    seen = set()
    for dir in path:
        normdir = os.path.normcase(dir)
        if not normdir in seen:
            seen.add(normdir)
            for thefile in files:
                name = os.path.join(dir, thefile)
                if _access_check(name, mode):
                    return name
    return None

def generate_input(structure):
    lines = ["P 1", "1"]
    lengths, angles = structure.lattice.lengths_and_angles
    cell = ["{:11.6f}".format(x) for x in lengths]
    cell.extend(["{:11.6f}".format(x) for x in angles])
    lines.append("".join(cell))
    for spec in structure:
        atom_line = ["{:8s}".format(spec.specie.symbol)]
        atom_line.extend(["{:10.6f}".format(x) for x in spec.frac_coords])
        lines.append("".join(atom_line))
    lines.extend(["end", "end"])
    return "\n".join(lines)


def parse_output(output):
    p = re.compile('(\D)+\d*')
    line = output.readline()
    while line:
        if line.startswith(" gamma/gamma(min) = 1.0000"):
            # This line has format: ' Cell parameters :   1.0802  1.0802  6.0100   89.974   89.974   86.583'
            tokens = [t for t in output.readline().split() if t]
            lattice = pymatgen.core.lattice.Lattice.from_parameters(*map(float, tokens[3:9]))

            species = list()
            coords = list()

            # These lines have format: ' A           0.00080  0.00080  0.66643                 1a'
            tokens = [t for t in output.readline().split() if t]
            while tokens:
                match = p.match(tokens[0])
                species.append(match.group(1))
                coords.append(map(float, tokens[1:4]))
                tokens = [t for t in output.readline().split() if t]

            return Structure(lattice, species, coords)

        line = output.readline()


def structure_tidy(structure):
    tidy = which("structure_tidy")
    if not tidy:
        return None
    tidy_exe_dir = os.path.dirname(os.path.abspath(tidy))
    os.environ["RIETAN"] = tidy_exe_dir

    with tempfile.NamedTemporaryFile('w') as f:
        f.write(generate_input(structure))
        f.flush()
        proc = subprocess.Popen([tidy, f.name], stdout=subprocess.PIPE)
        output = StringIO.StringIO()
        output.write(proc.communicate()[0])
        output.seek(0, 0)
        tidy_structure = parse_output(output)
        output.close()

    return tidy_structure