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
    gridCoords_.init(3, max(a, max(b, c)), 1, false);
    cout << "Dipoles" << endl;
    dipoles_.init(natoms, 7, 1, false);
    cout << "Grid contracted" << endl;
    gridContracted_.init(a*b*c, 4, 1, false);
}

void FieldMap::setupGrid(Frame *frame){
    /* create min and max initial values */
    gridBounds_(0, 0) = 1e6; gridBounds_(0, 1) = -1e6;
    gridBounds_(1, 0) = 1e6; gridBounds_(1, 1) = -1e6;
    gridBounds_(2, 0) = 1e6; gridBounds_(2, 1) = -1e6;
//    for(auto atom : frame->atoms_){
    for(int ii : frame->residues_[0].atoms){
        /* for each atom, compare min and max against coords */
            Atom *atom = &(frame->atoms_[ii]);
        gridBounds_(0, 0) = min(gridBounds_(0, 0), atom->coords[0]);
        gridBounds_(0, 1) = max(gridBounds_(0, 1), atom->coords[0]);
        gridBounds_(1, 0) = min(gridBounds_(1, 0), atom->coords[1]);
        gridBounds_(1, 1) = max(gridBounds_(1, 1), atom->coords[1]);
        gridBounds_(2, 0) = min(gridBounds_(2, 0), atom->coords[2]);
        gridBounds_(2, 1) = max(gridBounds_(2, 1), atom->coords[2]);
    }
    gridBounds_(0, 0) -= border_; gridBounds_(0, 1) += border_;
    gridBounds_(1, 0) -= border_; gridBounds_(1, 1) += border_;
    gridBounds_(2, 0) -= border_; gridBounds_(2, 1) += border_;
    /* set gridCentre */
    gridCentre_[0] = (gridBounds_(0, 1) - gridBounds_(0, 0)) / 2.;
    gridCentre_[1] = (gridBounds_(1, 1) - gridBounds_(1, 0)) / 2.;
    gridCentre_[2] = (gridBounds_(2, 1) - gridBounds_(2, 0)) / 2.;
    //cout << "Grid centre at: " << gridCentre_[0] << "," << gridCentre_[1] << "," << gridCentre_[2] << endl;
        gridCoords_.linspace(0, gridDims_[0], gridBounds_(0, 0), gridBounds_(0, 1));
        gridCoords_.linspace(1, gridDims_[1], gridBounds_(1, 0), gridBounds_(1, 1));
        gridCoords_.linspace(2, gridDims_[2], gridBounds_(2, 0), gridBounds_(2, 1));
//    for(int i=0; i<3; i++){
//        /* for x, y, z do linspace of grid coordinates */
//        gridCoords_.linspace(i, gridBounds_(i, 0), gridBounds_(i, 1));
//    }
//    cout << "Grid bounds" << endl;
//    cout << gridBounds_(0, 0) << " " << gridBounds_(0, 1) << endl;
//    cout << gridBounds_(1, 0) << " " << gridBounds_(1, 1) << endl;
//    cout << gridBounds_(2, 0) << " " << gridBounds_(2, 1) << endl;
//    cout << "Grid coords" << endl;
//    for(int i=0; i<gridDims_[0]; i++){
//        if(i%10 == 0) cout << endl;
//        cout << gridCoords_(0, i) << " ";
//    }
//    cout << endl;
//    for(int i=0; i<gridDims_[1]; i++){
//        if(i%10 == 0) cout << endl;
//        cout << gridCoords_(1, i) << " ";
//    }
//    cout << endl;
//    for(int i=0; i<gridDims_[2]; i++){
//        if(i%10 == 0) cout << endl;
//        cout << gridCoords_(2, i) << " ";
//    }
//    cout << endl;
//    cout << "Grid setup done" << endl;
}

void FieldMap::setupGridContracted(Frame *frame){
    /*RADMIN=50.0D0                                                     CHE01860        still inside x,y,z, loops
    187       DO 100 I=1,NATOMS                                                 CHE01870        for each atom
    188       VRAD = RADII(I)                                                   CHE01880
    189       DIST = (P1 - C(1,I))**2 + (P2 - C(2,I))**2 + (P3 - C(3,I))**2     CHE01890        calculate distance of grid point to atom
    190       DIST = DSQRT(DIST)                                                CHE01900
    191       IF (DIST .LT. VRAD) GOTO 210                                      CHE01910        if grid point is inside atom - skip
    192       IF (DIST .LT. RADMIN) RADMIN = DIST                               CHE01920        keep track of closest grid point
    193   100 CONTINUE                                                          CHE01930        do next atom
        194       IF (RADMIN .GT. RMAX) GOTO 210                                    CHE01940        if grid point is too far away from all atoms - skip
    195 C                                                                       CHE01950
    196 C        STORE POINTS (IN ATOMIC UNITS)                                 CHE01960
    197 C                                                                       CHE01970
    198       IPOINT = IPOINT + 1                                               CHE01980
    199       P(1,IPOINT) = P1                                                  CHE01990        store grid point if accepted
        200       P(2,IPOINT) = P2                                                  CHE02000
    201       P(3,IPOINT) = P3                                                  CHE02010
    202       IF (IPO(2) .EQ. 1)                                                CHE02020
    203      $ WRITE(IOUT,*) 'POINT ',IPOINT,' X,Y,Z ',P1,P2,P3                 CHE02030
    204   210 CONTINUE                                                          CHE02040
    205   200 CONTINUE
    */
    float radmin = 500.f;
    float dist = 0.f;
    float rmax = border_;    // use an rmax equal to border_ around molecule
    float vrad = 0.1f;       // reject if within 1A of an atom (inside atomic radius)
    float x, y, z;
    bool accepted;
    int accepted_count = 0, close_count = 0, far_count = 0;
    vector<float> coords(3);
    gridContracted_.appendedRows_ = 0;
#pragma omp parallel for
    for(int i=0; i < gridDims_[0]; i++){
        coords[0] = gridCoords_(0, i);
        for(int j=0; j < gridDims_[1]; j++){
            coords[1] = gridCoords_(1, j);
            for(int k=0; k < gridDims_[2]; k++){
                coords[2] = gridCoords_(2, k);
                radmin = 500.f;
                accepted = true;
                //for(Atom atom : frame->atoms_){
//                for(int ii : frame->residues_[0].atoms){
//                cout << i << j << k << endl;
                for(int ii=0; ii < frame->residues_[0].atoms.size(); ii++){
                    dist = sqrt(distSqr(frame->atoms_[ii].coords, coords[0], coords[1], coords[2]));
//                    cout << ii << endl;
//                    if(dist < border_) cout << dist << endl;
                    radmin = min(radmin, dist);
                    if(dist < vrad){
                        accepted = false;
//                        cout << "too close" << endl;
                        ++close_count;
                        break;
                    }
//                    if(dist < radmin) radmin = dist;
                }
                if(accepted && radmin > rmax){
//                    cout << "too far" << endl;
//                    cout << radmin << endl;
                    ++far_count;
                    continue;
                }
                if(accepted){
                    gridContracted_.append(coords);
                    ++accepted_count;
                }
            }
        }
    }
    numGridPoints_ = gridContracted_.appendedRows_;
    fieldMonopoleContracted_.reserve(numGridPoints_);
    cout << "Acc\tFar\tClose" << endl;
    cout << accepted_count << "\t" << far_count << "\t" << close_count << endl;
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

void FieldMap::calcFieldMonopolesContracted(Frame *frame){
    float inveps = 1. / (4 * M_PI * 8.854187817e-12);
    // inveps = 8.9875517873681e9
#pragma omp parallel for
    for(int i=0; i < numGridPoints_; i++) {
        for(Atom atom : frame->atoms_) {
            fieldMonopoleContracted_[i] += atom.charge /
                    distSqr(atom.coords, gridContracted_(i, 0), gridContracted_(i, 1), gridContracted_(i, 2));
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

