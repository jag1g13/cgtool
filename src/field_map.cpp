#include "field_map.h"

//#include <algorithm>
#include <cmath>
#include <iostream>

using std::min;
using std::max;
using std::vector;
using std::cout;
using std::endl;

FieldMap::FieldMap(){
}

FieldMap::FieldMap(const int a, const int b, const int c, const int natoms){
    gridDims_.reserve(3);
    gridDims_[0] = a; gridDims_[1] = b; gridDims_[2] = c;
    gridCentre_.reserve(3);
    cout << "Field Monopole" << endl;
    fieldMonopole_.init(a, b, c, false);
    cout << "Field Dipole" << endl;
    fieldDipole_.init(a, b, c, false);
    cout << "Grid bounds" << endl;
    gridBounds_.init(3, 2, 1, false);
    cout << "Coords" << endl;
    gridCoords_.init(max(a, max(b, c)), 3, 1, false);
    cout << "Dipoles" << endl;
    dipoles_.init(natoms, 7, 1, false);
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
    float x, y, z;
    #pragma omp parallel for
    for(int i=0; i < gridDims_[0]; i++){
        x = gridCoords_(i, 0);
        for(int j=0; j < gridDims_[1]; j++){
            y = gridCoords_(j, 1);
            for(int k=0; k < gridDims_[2]; k++){
                z = gridCoords_(k, 2);
                fieldMonopole_(i, j, k) = 0.;
                for(Atom atom : frame->atoms_){
                    fieldMonopole_(i, j, k) += atom.charge /
                            distSqr(atom.coords, x, y, z);
                }
            }
        }
    }
}

void FieldMap::calcFieldDipoles(Frame *frame) {
    float inveps = 1. / (4 * M_PI * 8.854187817e-12);
    //float inveps = 8.9875517873681e9;
    float x, y, z;
    float dip_angle = 0;
    #pragma omp parallel for
    for(int i=0; i < gridDims_[0]; i++){
        x = gridCoords_(i, 0);
        for(int j=0; j < gridDims_[1]; j++){
            y = gridCoords_(j, 1);
            for(int k=0; k < gridDims_[2]; k++){
                z = gridCoords_(k, 2);
                fieldDipole_(i, j, k) = 0.;
                int ii = 0;
                for(Atom atom : frame->atoms_){
                    fieldDipole_(i, j, k) += dipoles_(ii, 6) * cos(dip_angle) /
                            distSqr(atom.coords, x, y, z);
                    ii++;
                    //self.grid_dipole[i][j][k] += dipole[6] * cos(dip_angle)
                }
            }
        }
    }
}

/**
* \brief Return distance squared between two points.
*
* Originally used cmath pow - this version is much faster.
*/
float FieldMap::distSqr(float *coords, const float x, const float y, const float z) {
//    return (coords[0] - x)*(coords[0] - x) +
//            (coords[1] - y)*(coords[1] - y) +
//            (coords[2] - z)*(coords[2] - z);
    float tmpx = coords[0] - x;
    float tmpy = coords[1] - y;
    float tmpz = coords[2] - z;
    return tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
}

