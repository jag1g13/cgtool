#!/usr/bin/env python
__author__ = 'James Graham'

import struct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import axes3d



def read_file(filename="/home/james/cgtool/build/field200_cg.csv"):
    with open(filename) as infile:
        points = []
        xmin, ymin, zmin = 1000, 1000, 1000
        xmax, ymax, zmax = -999, -999, -999
        for line in infile:
            # print(len(line))
            # print(line)
            (x, y, z, f) = struct.unpack("12s12s12s12s", line[0:-1])
            (x, y, z, f) = map(float, (x, y, z, f))
            points.append((x, y, z, f))
            xmin = min(x, xmin)
            xmax = max(x, xmax)
            ymin = min(y, ymin)
            ymax = max(y, ymax)
            zmin = min(z, zmin)
            zmax = max(z, zmax)

        grid = 50
        field_map = np.zeros([grid, grid, grid])

        xstep = (xmax - xmin) / (grid-1)
        ystep = (ymax - ymin) / (grid-1)
        zstep = (zmax - zmin) / (grid-1)
        for (x, y, z, f) in points:
            x_bin = (x - xmin) / xstep
            y_bin = (y - ymin) / ystep
            z_bin = (z - zmin) / zstep
            # int loc = int((val - min_) / step_);
            field_map[x_bin, y_bin, z_bin] += f

        # print(field_map)
        xs = np.linspace(xmin, xmax, grid)
        ys = np.linspace(ymin, ymax, grid)
        XS, YS = np.meshgrid(xs, ys)
        fig = plt.figure()
        im = plt.imshow(field_map[:, :, grid/2], interpolation='bilinear', origin='lower', cmap=cm.gray, extent=(xmin, xmax, ymin, ymax))
        # CS = plt.contour(XS, YS, field_map[:, :, grid/2])
        # ax = fig.add_subplot(111, projection='3d')
        # ax.plot_wireframe(XS, YS, field_map[:, :, grid/2])


def main():
    read_file("/home/james/projects/cgtool/build/field200_aa.csv")
    read_file("/home/james/projects/cgtool/build/field200_cg.csv")
    plt.show()


if __name__ == "__main__":
    main()
