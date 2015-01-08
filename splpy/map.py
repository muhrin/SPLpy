"""
Module that processes output from stools' smap.
"""

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014, Martin Uhrin"
__version__ = "0.0.1"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "Dec 30, 2014"

import math
import os
import re

import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches

from shapely.geometry import Polygon

OUTPUTTERS = ['matplotlib', 'latex']

def parse_point(point):
    nums = point[1:-1].split(',')
    return float(nums[0]), float(nums[1])


class MatplotlibOutputter:
    def out(self, map_file, output_name=None, settings=None):
        fig, ax = plt.subplots()
        output_file = output_name if output_name else "map"
        output_file += ".pdf"

        line = map_file.readline()
        while line:
            line = line.rstrip(os.linesep)
            if line != '':
                self.draw_face(ax, map_file, line, settings)
            line = map_file.readline()

        ax.set_aspect('equal', 'box')
        # Set the x and y range to be exactly the data limits
        ax.set_xlim([ax.dataLim.xmin, ax.dataLim.xmax])
        ax.set_ylim([ax.dataLim.ymin, ax.dataLim.ymax])

        if settings:
            if "title" in settings:
                plt.title(settings["title"], fontsize=18)
            if "xlabel" in settings:
                plt.xlabel(settings["xlabel"], fontsize=17)
            if "ylabel" in settings:
                plt.ylabel(settings["ylabel"], fontsize=17)

        plt.draw()
        plt.savefig(output_file, bbox_inches='tight')
        plt.close()

    def draw_face(self, axes, map_file, label, settings):
        line = map_file.readline().rstrip(os.linesep)
        if line == 'Points':
            self.draw_points(map_file, label, settings)
            self.draw_path(axes, map_file, label, settings)

    def draw_points(self, map_file, label, settings):

        line = map_file.readline().rstrip(os.linesep)
        x = list()
        y = list()
        while line:
            if line == 'Path':
                break
            pt = parse_point(line)
            x.append(pt[0])
            y.append(pt[1])

            line = map_file.readline().rstrip(os.linesep)

        plt.scatter(x, y, color=self.get_property(label, 'color', settings), s=0.3, zorder=1)

    def draw_path(self, axes, map_file, label, settings):
        line = map_file.readline().rstrip(os.linesep)
        codes = list()
        coords = list()
        poly_coords = list()
        while line and len(line) >= 6:
            command = line[:6]
            if command == 'MOVETO':
                codes.append(mpath.Path.MOVETO)
                pt = parse_point(line[7:])
                coords.append(pt)
                poly_coords.append(pt)
            elif command == 'LINETO':
                codes.append(mpath.Path.LINETO)
                pt = parse_point(line[7:])
                coords.append(pt)
                poly_coords.append(pt)
            elif command == 'CURVE4':
                codes.append(mpath.Path.CURVE4)
                codes.append(mpath.Path.CURVE4)
                codes.append(mpath.Path.CURVE4)

                match = re.search(r'(\(.+\)), (\(.+\)), (\(.+\))', line[7:])
                coords.append(parse_point(match.group(1)))
                coords.append(parse_point(match.group(2)))
                end_pt = parse_point(match.group(3))
                coords.append(end_pt)
                poly_coords.append(end_pt)
            else:
                return
            line = map_file.readline().rstrip(os.linesep)

        face_path = mpath.Path(coords, codes, closed=True)
        face_patch = mpatches.PathPatch(face_path, fill=False, color='black', linewidth=0.5)
        axes.add_patch(face_patch)

        poly = Polygon(poly_coords)
        if poly.is_valid and poly.area > 0.09:
            rep_pt = poly.representative_point()
            fontsize = min(28 * math.sqrt(poly.area) + 2.0, 22)
            color = self.get_property(label, 'color', settings)
            plt.text(rep_pt.x, rep_pt.y, label, size=fontsize,
                     horizontalalignment='center', verticalalignment='center',
                     bbox=dict(facecolor=color, edgecolor=color, boxstyle='round', alpha=0.75))

    def get_property(self, label, prop, settings):
        if not settings or 'labels' not in settings:
            return None

        all_properties = settings['labels']
        if label in all_properties:
            properties = all_properties[label]
            if prop in properties:
                return properties[prop]
        return 'gray'


class LatexOutputter:
    def out(self, map_file, settings):
        line = map_file.readline().rstrip(os.linesep)
        while line is not None:
            if line != '':
                self.draw_face(map_file, line)
            line = map_file.readline().rstrip(os.linesep)

    @classmethod
    def draw_points(cls, map_file, label):
        print('\\addplot [only marks] table {')

        finished_points = False
        line = map_file.readline().rstrip(os.linesep)
        while line:
            if line == 'Path':
                finished_points = True
                break
            x, y = parse_point(line)
            print('{} {}'.format(x, y))
            line = map_file.readline().rstrip(os.linesep)

        print('};')
        return finished_points

    @classmethod
    def draw_path(cls, map_file, label):
        print('\draw '),
        line = map_file.readline().rstrip(os.linesep)
        while line:
            command = line[:6]
            if command == 'MOVETO':
                print('(axis cs: {})'.format(line[8:-1])),
            elif command == 'LINETO':
                print('-- (axis cs: {})'.format(line[8:-1])),
            elif command == 'CURVE4':
                match = re.search(r'\((.+)\), \((.+)\), \((.+)\)', line[7:])
                print('.. controls (axis cs: {}) and (axis cs: {}) .. (axis cs: {})').format(match.group(1),
                                                                                             match.group(2),
                                                                                             match.group(3)),
            else:
                print(';')
                return
            line = map_file.readline().rstrip(os.linesep)

        print(';')

    @classmethod
    def draw_face(cls, map_file, label):
        line = map_file.readline().rstrip(os.linesep)
        if line == 'Points':
            if cls.draw_points(map_file, label):
                # Now we should be ready to draw the path
                cls.draw_path(map_file, label)


def draw_map(map_file, outputter, settings=None, output_name=None):
    if outputter == 'matplotlib':
        out = MatplotlibOutputter()
    elif outputter == 'latex':
        out = LatexOutputter()
    else:
        raise ValueError("Unknown map outputter: {}".format(outputter))

    out.out(map_file, output_name, settings)
