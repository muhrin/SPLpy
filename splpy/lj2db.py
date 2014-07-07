#!/usr/bin/python

__author__ = 'Martin Uhrin'

import argparse

from pymatgen.apps.borg.queen import BorgQueen
from splpy.db.lj_assimilate import LjToDbTaskDrone

parser = argparse.ArgumentParser(
    description='Insert Lennard-Jones structures into a database with potential parameters')
parser.add_argument('assimilate_path', default='.', nargs='?',
                    help='The path to the potparams file or a directory with multiple potparams files')


args = parser.parse_args()
drone = LjToDbTaskDrone()
queen = BorgQueen(drone, args.assimilate_path, 2)
#entries = queen.get_data()
