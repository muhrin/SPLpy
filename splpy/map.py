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
import StringIO
import subprocess
import tempfile

import matplotlib
# Make sure matplotlib doesn't try to use an xserver requiring backend by default
# Must be before importing matplotlib.pyplot or pylab!
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches

from shapely.geometry import Polygon
from shapely.geometry import point

OUTPUTTERS = ['matplotlib', 'latex']


def parse_point(point):
    nums = point[1:-1].split(',')
    return float(nums[0]), float(nums[1])


class MatplotlibOutputter:
    def out(self, map_file, output_name=None, options=None, mask=None):
        settings = self._init_settings(options)

        fig, ax = plt.subplots()
        output_file = output_name if output_name else "map"
        output_file += ".pdf"

        self._draw(map_file, ax, settings)
        if mask:
            self.draw_mask(mask, ax, settings)

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

        plt.rc('pdf', fonttype=42)
        plt.rc('ps', fonttype=42)
        plt.draw()
        plt.savefig(output_file, bbox_inches='tight', dpi=600)
        plt.close()

    def draw_face(self, axes, map_file, label, settings):
        line = map_file.readline().rstrip(os.linesep)
        if line == 'Points':
            x, y = self._parse_points(map_file)
            self.draw_points(x, y, self.get_property(label, 'color', settings), 2, settings)

        codes, coords = self.parse_path(map_file)
        self._draw_path(axes, coords, codes, settings, label=label)

    def draw_points(self, x, y, color, size, settings):
        path = plt.scatter(x, y, color=color, s=size, lw=0, zorder=1)
        path.set_rasterized(settings['rasterized'])


    def draw_mask(self, map_file, axes, settings):
        found_mask = False
        for line in map_file:
            line = line.rstrip(os.linesep)

            if line == 'False':
                found_mask = True
            if found_mask and line == 'Path':
                codes, coords = self.parse_path(map_file)
                face_path = mpath.Path(coords, codes, closed=True)
                mask_settings = settings['mask']
                face_patch = mpatches.PathPatch(face_path, fill=True, facecolor=mask_settings['color'],
                                                edgecolor='#000000', linewidth=0,
                                                alpha=mask_settings['alpha'], zorder=5)
                axes.add_patch(face_patch)
                found_mask = False

    def _parse_points(self, map_file):
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
        return x, y

    def parse_path(self, map_file):
        codes = list()
        coords = list()
        line = map_file.readline().rstrip(os.linesep)
        while line and len(line) >= 6:
            command = line[:6]
            if command == 'MOVETO':
                codes.append(mpath.Path.MOVETO)
                pt = parse_point(line[7:])
                coords.append(pt)
            elif command == 'LINETO':
                codes.append(mpath.Path.LINETO)
                pt = parse_point(line[7:])
                coords.append(pt)
            elif command == 'CURVE4':
                codes.append(mpath.Path.CURVE4)
                codes.append(mpath.Path.CURVE4)
                codes.append(mpath.Path.CURVE4)

                match = re.search(r'(\(.+\)), (\(.+\)), (\(.+\))', line[7:])
                for i in range(1, 4):
                    coords.append(parse_point(match.group(i)))
            else:
                return
            line = map_file.readline().rstrip(os.linesep)
        return codes, coords

    def _draw_path(self, axes, coords, codes, settings, fill=False, color='black', linewidth=0.5, label=None):
        face_path = mpath.Path(coords, codes, closed=True)
        face_patch = mpatches.PathPatch(face_path, fill=fill, color=color, linewidth=linewidth)
        axes.add_patch(face_patch)

        if label:
            poly_coords = list()
            num_curve_codes = 0
            for pt, code in zip(coords, codes):
                if code is 'CURVE4':
                    num_curve_codes += 1
                    if num_curve_codes == 3:
                        poly_coords.append(pt)
                        num_curve_codes = 0
                else:
                    poly_coords.append(pt)
                    num_curve_codes = 0

            rep_pt = None
            poly = Polygon(poly_coords)
            num_points = len(poly_coords)
            if poly.is_valid and poly.area > 0.075:
                rep_pt = poly.representative_point()
                fontsize = min(26 * math.sqrt(poly.area) + 2.0, 20)
            elif num_points > 45:
                sum_x = 0.0
                sum_y = 0.0
                for coord in coords:
                    sum_x += coord[0]
                    sum_y += coord[1]
                rep_pt = point.Point([sum_x / float(num_points), sum_y / float(num_points)])
                fontsize = min(float(num_points) / 12.0, 20)
            if rep_pt:
                color = self.get_property(label, 'color', settings)
                plt.text(rep_pt.x, rep_pt.y, self.get_label_string(label, settings), size=fontsize,
                         horizontalalignment='center', verticalalignment='center',
                         bbox=dict(facecolor=color, edgecolor=color, boxstyle='round', alpha=0.75))

    def get_property(self, label, prop, settings):
        value = None

        if settings and 'labels' in settings:
            all_properties = settings['labels']
            if label in all_properties:
                properties = all_properties[label]
                if prop in properties:
                    value = properties[prop]

        # Defaults
        if not value:
            if prop == 'color':
                value = 'gray'

        return value

    def get_label_string(self, label, settings):
        label_str = self.get_property(label, 'latex_label', settings)

        # Default
        if not label_str:
            label_str = label

        return label_str

    def _draw(self, map_file, ax, settings=None):
        for line in map_file:
            line = line.rstrip(os.linesep)
            if line != '':
                self.draw_face(ax, map_file, line, settings)

    def _init_settings(self, settings):
        # Put the defaults here
        defaults = {
            'razterized': False,
            'mask': {
                'color': 'gray',
                'alpha': 0.4,
            }
        }
        defaults.update(settings)
        return defaults


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


def generate_map(map_points):
    output = StringIO.StringIO()
    with tempfile.NamedTemporaryFile('w') as f:
        for pt in map_points:
            f.write('{},{},{}\n'.format(pt[0], pt[1], pt[2]))
        f.flush()

        proc = subprocess.Popen(["smap", f.name], stdout=subprocess.PIPE)
        output.write(proc.communicate()[0])

    output.seek(0, 0)
    return output


def draw_map(map_file, outputter, settings=None, output_name=None, mask=None):
    if outputter == 'matplotlib':
        out = MatplotlibOutputter()
    elif outputter == 'latex':
        out = LatexOutputter()
    else:
        raise ValueError("Unknown map outputter: {}".format(outputter))

    out.out(map_file, output_name, settings, mask)
