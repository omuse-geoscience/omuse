#!/usr/bin/env python


from omuse.community.adcirc.read_grid import adcirc_grid_reader
import numpy
from matplotlib import pyplot

class adcirc_grid_viewer():

    def __init__(self, filename, coordinates):
        r = self.r = adcirc_grid_reader(filename, coordinates)
        r.read_grid()

        # r.p contains the lines of node definitions, strip last column
        self.x = r.p[:,0]
        self.y = r.p[:,1]

        # r.t contains the connectivity of nodes to form triangles, strip fist column
        self.triangles = r.t[:,1:] -1 

    def plot(self):
        pyplot.triplot(self.x, self.y, self.triangles, 'k-')
        pyplot.show()


if __name__ == "__main__":
    import sys
    import os.path

    usage_str = "Usage: view_adcirc_grid.py filename coordinates=[spherical|cartesian]"

    total = len(sys.argv)
    if len(sys.argv) != 3:
        sys.exit(usage_str)

    filename = sys.argv[1]
    if not os.path.isfile(filename):
        sys.exit("Error: No such file " + filename)

    coordinates = sys.argv[2]
    if not (coordinates == 'spherical' or coordinates == 'cartesian'):
        sys.exit(usage_str)

    v = adcirc_grid_viewer(filename, coordinates)

    v.plot()

    raw_input()


