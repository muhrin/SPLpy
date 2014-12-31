
import argparse
import os
import yaml

import splpy.map

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Draw a map from smap raw output')
    parser.add_argument("map", help='The smap file to draw')
    parser.add_argument("--settings", dest="settings", type=str, default=None, required=False,
                        help="The yaml settings file with information about the plot")

    args = parser.parse_args()

    if args.settings and os.path.exists(args.settings):
        with open(args.settings, 'r') as settings_file:
            settings = yaml.load(settings_file)

    with open(args.map, 'r') as map_file:
        splpy.map.draw_map(map_file, 'matplotlib', settings, os.path.splitext(args.map)[0])

