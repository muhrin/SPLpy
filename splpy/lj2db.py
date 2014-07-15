#!/usr/bin/python

__author__ = 'Martin Uhrin'

import argparse

from pymatgen.apps.borg.queen import BorgQueen
from splpy.lj_assimilate import LjToDbTaskDrone

parser = argparse.ArgumentParser(
    description='Insert Lennard-Jones structures into a database with potential parameters')
parser.add_argument('assimilate_path', default='.', nargs='?',
                    help='The folder containing potparams files')
parser.add_argument('-host', default="127.0.0.1", help='The MongoDB host to insert to')
parser.add_argument('-tags', help='Tag this set of structures with one or more tags, separate with commas.')

args = parser.parse_args()

tags = args.tags.split(',')

drone = LjToDbTaskDrone(args.host, tags=tags)
queen = BorgQueen(drone, args.assimilate_path, 2)
#entries = queen.get_data()
