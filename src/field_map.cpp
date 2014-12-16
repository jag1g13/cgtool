#include <stdexcept>
#include <algorithm>
#include <iostream>

#include "arrays.h"
#include "field_map.h"

using std::min;
using std::max;
using std::vector;
using std::cout;
using std::endl;

FieldMap::FieldMap(){
}

FieldMap::FieldMap(const int a, const int b, const int c){
    gridDims_.reserve(3);
    gridDims_[0] = a; gridDims_[1] = b; gridDims_[2] = c;
//    fieldMonopole_ = alloc_float_3d(a, b, c);
//    if(fieldMonopole_ == NULL) throw std::runtime_error("Array alloc error");
    //fieldMonopole_ = ArrayFloat(gridDims_);
    cout << "Field Monopole" << endl;
    //fieldMonopole_ = ArrayFloat(a, b, c);
    fieldMonopole_.init(a, b, c);
//    fieldDipole_ = alloc_float_3d(a, b, c);
    //fieldDipole_ = ArrayFloat(gridDims_);
    cout << "Field Dipole" << endl;
    //fieldDipole_ = ArrayFloat(a, b, c);
    fieldDipole_.init(a, b, c);
//    if(fieldDipole_ == NULL) throw std::runtime_error("Array alloc error");
//    gridBounds_ = alloc_float_2d(3, 2);
    //vector<int> tmp{3, 2};
    //gridBounds_ = ArrayFloat(tmp);
    cout << "Grid bounds" << endl;
    //gridBounds_ = ArrayFloat(3, 2, 1);
    gridBounds_.init(3, 2, 1);
//    if(gridBounds_ == NULL) throw std::runtime_error("Array alloc error");
//    gridCoords_ = alloc_float_2d(3, max(a, max(b, c)));
    //tmp = {3, max(a, max(b, c))};
    //gridCoords_ = ArrayFloat(tmp);
    cout << "Coords" << endl;
    //gridCoords_ = ArrayFloat(3, max(a, max(b, c)), 1);
    gridCoords_.init(3, max(a, max(b, c)), 1);
//    if(gridCoords_ == NULL) throw std::runtime_error("Array alloc error");
}

void FieldMap::setupGrid(Frame *frame){
    /* create min and max initial values */
//    gridBounds_[0][0] = 1e6; gridBounds_[0][1] = -1e6;
//    gridBounds_[1][0] = 1e6; gridBounds_[1][1] = -1e6;
//    gridBounds_[2][0] = 1e6; gridBounds_[2][1] = -1e6;
    gridBounds_(0, 0) = 1e6; gridBounds_(0, 1) = -1e6;
    gridBounds_(1, 0) = 1e6; gridBounds_(1, 1) = -1e6;
    gridBounds_(2, 0) = 1e6; gridBounds_(2, 1) = -1e6;
    for(auto atom : frame->atoms_){
        /* for each atom, compare min and max against coords */
//        gridBounds_[0][0] = min(gridBounds_[0][0], atom.coords[0]);
//        gridBounds_[0][1] = max(gridBounds_[0][1], atom.coords[0]);
//        gridBounds_[1][0] = min(gridBounds_[1][0], atom.coords[1]);
//        gridBounds_[1][1] = max(gridBounds_[1][1], atom.coords[1]);
//        gridBounds_[2][0] = min(gridBounds_[2][0], atom.coords[2]);
//        gridBounds_[2][1] = max(gridBounds_[2][1], atom.coords[2]);
        gridBounds_(0, 0) = min(gridBounds_(0, 0), atom.coords[0]);
        gridBounds_(0, 1) = max(gridBounds_(0, 1), atom.coords[0]);
        gridBounds_(1, 0) = min(gridBounds_(1, 0), atom.coords[1]);
        gridBounds_(1, 1) = max(gridBounds_(1, 1), atom.coords[1]);
        gridBounds_(2, 0) = min(gridBounds_(2, 0), atom.coords[2]);
        gridBounds_(2, 1) = max(gridBounds_(2, 1), atom.coords[2]);
    }
//    gridBounds_[0][0] -= border_; gridBounds_[0][1] += border_;
//    gridBounds_[1][0] -= border_; gridBounds_[1][1] += border_;
//    gridBounds_[2][0] -= border_; gridBounds_[2][1] += border_;
    gridBounds_(0, 0) -= border_; gridBounds_(0, 1) += border_;
    gridBounds_(1, 0) -= border_; gridBounds_(1, 1) += border_;
    gridBounds_(2, 0) -= border_; gridBounds_(2, 1) += border_;
    for(int i=0; i<3; i++){
        /* for x, y, z do linspace of grid coordinates */
//        linspace_1d(gridCoords_[0], gridBounds_[0][0], gridBounds_[0][1], gridDims_[0]);
//        linspace_1d(gridCoords_[1], gridBounds_[1][0], gridBounds_[1][1], gridDims_[1]);
//        linspace_1d(gridCoords_[2], gridBounds_[2][0], gridBounds_[2][1], gridDims_[2]);
    }
}

/*
import numpy as np
import matplotlib.pyplot as plt
# from scipy import optimize
from math import sqrt
from math import cos

# from frame import Frame


class FieldMap:
    """
    Calclates, stores and manipulates the electric field around a molecule, given from a Frame instance
    """
    def __init__(self, frame):
        self.border = 1       # disance from molecule that the grid extends in nm
        # self.grid_dim = [3, 75, 75]       # looks okay as an image
        # self.grid_dim = [3, 125, 125]     # looks good as an image
        self.grid_dim = [5, 5, 5]     # looks bad as an image
        self.grid_monopole = np.zeros(self.grid_dim)
        self.grid_dipole = np.zeros(self.grid_dim)
        self.dipoles = np.ndarray([len(frame.atoms), 7])  # for each atom store dipole coords, vector and magnitude
        print("Grid:", self.grid_dim)
        print(str(self.grid_dim[0] * self.grid_dim[1] * self.grid_dim[2]) + " grid points")

    def setup_grid(self, frame):
        """
        Need to check the bounds of the grid with every new frame
        :param frame: Frame object containing the xtc frame to read
        :return: Nothing
        """
        xmax, ymax, zmax = float("-inf"), float("-inf"), float("-inf")
        xmin, ymin, zmin = float("inf"), float("inf"), float("inf")
        for atom in frame.atoms:
            xmax = max(xmax, atom.loc[0])
            xmin = min(xmin, atom.loc[0])
            ymax = max(ymax, atom.loc[1])
            ymin = min(ymin, atom.loc[1])
            zmax = max(zmax, atom.loc[2])
            zmin = min(zmin, atom.loc[2])
        xmax += self.border
        ymax += self.border
        zmax += self.border
        xmin -= self.border
        ymin -= self.border
        zmin -= self.border
        self.grid_x = np.linspace(xmin, xmax, self.grid_dim[0])
        self.grid_y = np.linspace(ymin, ymax, self.grid_dim[1])
        self.grid_z = np.linspace(zmin, zmax, self.grid_dim[2])

    def calc_field_monopoles(self, frame):
        """
        Calculate the electric field over the grid from atomic charges
        :param frame:
        :return: Nothing, store result in self.grid_monopole
        """
        inveps = 1. / (4 * np.pi * 8.854187817e-12)     # I don't multiply by this, would just cancel out anyway
        # inveps = 8.9875517873681e9
        for i in xrange(self.grid_dim[0]):
            for j in xrange(self.grid_dim[1]):
                for k in xrange(self.grid_dim[2]):
                    self.grid_monopole[i][j][k] = 0.
                    for atom in frame.atoms:
                        self.grid_monopole[i][j][k] += atom.charge /\
                                                       self.dist_sqr(atom.loc[0], atom.loc[1], atom.loc[2],
                                                                     i, j, k)
                    # self.grid[i][j][k] *= inveps

    def calc_field_dipoles(self):
        """
        Calculate the electric field over the grid from point dipoles
        :return: Nothing, store result in self.grid_dipole
        """
        inveps = 1. / (4 * np.pi * 8.854187817e-12)     # I don't multiply by this, would just cancel out anyway
        # inveps = 8.9875517873681e9
        for i in xrange(self.grid_dim[0]):
            for j in xrange(self.grid_dim[1]):
                for k in xrange(self.grid_dim[2]):
                    self.grid_dipole[i][j][k] = 0.
                    for dipole in self.dipoles:
                        #self.grid_monopole[i][j][k] += atom.charge / self.dist_sqr(atom, i, j, k)
                        self.grid_dipole[i][j][k] += dipole[6] * cos(dip_angle) /\
                                                     dist_sqr(dipole[0], dipole[1], dipole[2],
                                                              i, j, k)
                    # self.grid[i][j][k] *= inveps

    def dist_sqr(self, x, y, z, i, j, k):
        return (x-self.grid_x[i])**2 + \
               (y-self.grid_y[j])**2 + \
               (z-self.grid_z[k])**2

    def plot(self, x):
        # plt.contourf(self.grid_y, self.grid_z, self.grid[x])
        plt.pcolor(self.grid_y, self.grid_z, self.grid_monopole[x], vmin=-1, vmax=1)
        # plt.imshow(self.grid[x], extent=[self.grid_y[0], self.grid_y[-1], self.grid_z[0], self.grid_z[-1]])
        plt.show()

    def __repr__(self):
        # print("Grid:", self.grid_dim)
        # print(self.grid)
        return "Grid: " + repr(self.grid_dim) + "\n" +\
               str(self.grid_dim[0] * self.grid_dim[1] * self.grid_dim[2]) + " grid points\n" +\
               repr(self.grid_monopole)
*/
