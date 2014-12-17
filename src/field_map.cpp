#include "field_map.h"

#include <algorithm>
#include <iostream>

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
    gridCentre_.reserve(3);
    cout << "Field Monopole" << endl;
    fieldMonopole_.init(a, b, c, true);
    cout << "Field Dipole" << endl;
    fieldDipole_.init(a, b, c, true);
    cout << "Grid bounds" << endl;
    gridBounds_.init(3, 2, 1, true);
    cout << "Coords" << endl;
    gridCoords_.init(3, max(a, max(b, c)), 1, true);
}

void FieldMap::setupGrid(Frame *frame){
    /* create min and max initial values */
    gridBounds_(0, 0) = 1e6; gridBounds_(0, 1) = -1e6;
    gridBounds_(1, 0) = 1e6; gridBounds_(1, 1) = -1e6;
    gridBounds_(2, 0) = 1e6; gridBounds_(2, 1) = -1e6;
    for(auto atom : frame->atoms_){
        /* for each atom, compare min and max against coords */
        gridBounds_(0, 0) = min(gridBounds_(0, 0), atom.coords[0]);
        gridBounds_(0, 1) = max(gridBounds_(0, 1), atom.coords[0]);
        gridBounds_(1, 0) = min(gridBounds_(1, 0), atom.coords[1]);
        gridBounds_(1, 1) = max(gridBounds_(1, 1), atom.coords[1]);
        gridBounds_(2, 0) = min(gridBounds_(2, 0), atom.coords[2]);
        gridBounds_(2, 1) = max(gridBounds_(2, 1), atom.coords[2]);
    }
    gridBounds_(0, 0) -= border_; gridBounds_(0, 1) += border_;
    gridBounds_(1, 0) -= border_; gridBounds_(1, 1) += border_;
    gridBounds_(2, 0) -= border_; gridBounds_(2, 1) += border_;
    /* set gridCentre */
    gridCentre_[0] = (gridBounds_(0, 1) - gridBounds_(0, 0)) / 2.;
    gridCentre_[1] = (gridBounds_(1, 1) - gridBounds_(1, 0)) / 2.;
    gridCentre_[2] = (gridBounds_(2, 1) - gridBounds_(2, 0)) / 2.;
    //cout << "Grid centre at: " << gridCentre_[0] << "," << gridCentre_[1] << "," << gridCentre_[2] << endl;
    for(int i=0; i<3; i++){
        /* for x, y, z do linspace of grid coordinates */
        gridCoords_.linspace(0, gridBounds_(0, 0), gridBounds_(0, 1));
        gridCoords_.linspace(1, gridBounds_(1, 0), gridBounds_(1, 1));
        gridCoords_.linspace(2, gridBounds_(2, 0), gridBounds_(2, 1));
    }
}

void FieldMap::calcFieldMonopoles(Frame *frame){
    float inveps = 1. / (4 * M_PI * 8.854187817e-12);
    // inveps = 8.9875517873681e9
    for(int i=0; i < gridDims_[0]; i++){
        for(int j=0; j < gridDims_[1]; j++){
            for(int k=0; k < gridDims_[2]; k++){
                fieldMonopole_(i, j, k) = 0.;
                for(Atom atom : frame->atoms_){
                    fieldMonopole_(i, j, k) += atom.charge /
                            distSqr(atom.coords, i, j, k);
                }
            }
        }
    }
}

float FieldMap::distSqr(float *coords, int i, int j, int k) {
    return pow((coords[0] - gridCoords_(0, i)), 2.f) +
            pow((coords[1] - gridCoords_(1, j)), 2.f) +
            pow((coords[2] - gridCoords_(2, k)), 2.f);
}

/*
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
*/
