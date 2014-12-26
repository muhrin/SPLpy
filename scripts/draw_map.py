import argparse
import os
import re
import yaml

import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches


def parse_point(point):
    nums = point[1:-1].split(',')
    return float(nums[0]), float(nums[1])


class MatplotlibOutputter:
    def __init__(self, settings):
        self.settings = settings

    def out(self, map_file):
        fig, ax = plt.subplots()

        line = map_file.readline()
        while line:
            line = line.rstrip(os.linesep)
            if line != '':
                self.draw_face(ax, map_file, line)
            line = map_file.readline()
        plt.legend()
        plt.show()

    def draw_face(self, axes, map_file, label):
        line = map_file.readline().rstrip(os.linesep)
        if line == 'Points':
            x, y = self.get_points(map_file, label)
            # Now we should be ready to draw the path
            self.draw_path(axes, map_file, label)
            # Finally draw the points on top of the path
            # plt.scatter(x, y, color=self.get_property(label, 'color'), s=0.5)
            plt.scatter(x, y, color='black', s=0.5)

    def get_points(self, map_file, label):

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

    def draw_path(self, axes, map_file, label):
        line = map_file.readline().rstrip(os.linesep)
        codes = list()
        coords = list()
        while line and len(line) >= 6:
            command = line[:6]
            if command == 'MOVETO':
                codes.append(mpath.Path.MOVETO)
                coords.append(parse_point(line[7:]))
            elif command == 'LINETO':
                codes.append(mpath.Path.LINETO)
                coords.append(parse_point(line[7:]))
            elif command == 'CURVE4':
                codes.append(mpath.Path.CURVE4)
                codes.append(mpath.Path.CURVE4)
                codes.append(mpath.Path.CURVE4)

                match = re.search(r'(\(.+\)), (\(.+\)), (\(.+\))', line[7:])
                coords.append(parse_point(match.group(1)))
                coords.append(parse_point(match.group(2)))
                coords.append(parse_point(match.group(3)))
            else:
                return
            line = map_file.readline().rstrip(os.linesep)

        face_path = mpath.Path(coords, codes, closed=True)
        face_patch = mpatches.PathPatch(face_path, alpha=0.5, facecolor='white', color=self.get_property(label, 'color'), label=label)
        axes.add_patch(face_patch)

    def get_property(self, label, prop):
        if not self.settings or 'labels' not in self.settings:
            return None

        all_properties = self.settings['labels']
        if label in all_properties:
            properties = all_properties[label]
            if prop in properties:
                return properties[prop]
        return None


class LatexOutputter:
    def out(self, map_file):
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Draw a map from smap raw output')
    parser.add_argument("map", help='The smap file to draw')
    parser.add_argument("--settings", dest="settings", type=str, default=None, required=False,
                        help="The yaml settings file with information about the plot")

    args = parser.parse_args()

    if args.settings and os.path.exists(args.settings):
        with open(args.settings, 'r') as settings_file:
            settings = yaml.load(settings_file)

    outputter = MatplotlibOutputter(settings)
    with open(args.map, 'r') as map_file:
        outputter.out(map_file)

