//
// Created by james on 28/04/15.
//

#include "membrane.h"

#include <cmath>
#include <iostream>

#include "small_functions.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::map;

Membrane::Membrane(const vector<Residue> &residues){
    residues_ = residues;
}

void Membrane::sortBilayer(const Frame &frame, const int blocks){
//TODO fix this so that ref_atom can be used by name
//    const int ref_atom = frame.nameToNum_.at(residue_.ref_atom);

    // Copy box from Frame - assume orthorhombic
    box_[0] = frame.box_[0][0];
    box_[1] = frame.box_[1][1];
    box_[2] = frame.box_[2][2];

    // Find middle of membrane in z coord
    // Do this in blocks to account for curvature - more curvature needs more blocks
    Array block_avg_z(blocks, blocks);
    Array block_tot_residues(blocks, blocks);

    for(const Residue &res : residues_){
        if(res.ref_atom < 0) continue;
        numLipids_ += res.num_residues;
        for(int i = 0; i < res.num_residues; i++){
            const int num = res.ref_atom + i * res.num_atoms + res.start;
            const int x = int(fmax((frame.x_[num][0] * blocks / box_[0]), 0.));
            const int y = int(fmin((frame.x_[num][1] * blocks / box_[1]), blocks-1));
            block_avg_z(x, y) += frame.x_[num][2];
            block_tot_residues(x, y)++;
        }
    }
    block_avg_z.elementDivide(block_tot_residues);

    // Separate membrane into upper and lower
    for(const Residue &res : residues_){
        if(res.ref_atom < 0) continue;
        int num_in_leaflet[2] = {0, 0};
        for(int i = 0; i < res.num_residues; i++){
            const int num = res.ref_atom + i * res.num_atoms + res.start;
            const int x = int(frame.x_[num][0] * blocks / box_[0]);
            const int y = int(frame.x_[num][1] * blocks / box_[1]);
            const double z = frame.x_[num][2];

            if(z < block_avg_z(x, y)){
                lowerHeads_.push_back(num);
                num_in_leaflet[0]++;
            } else{
                upperHeads_.push_back(num);
                num_in_leaflet[1]++;
            }
        }
        printf("%5s: %'4d lower, %'4d upper\n",
               res.resname.c_str(), num_in_leaflet[0], num_in_leaflet[1]);
    }
}

void Membrane::thickness(const Frame &frame, const bool with_reset){
    if(with_reset) reset();

    // Copy box from Frame - assume orthorhombic
    box_[0] = frame.box_[0][0];
    box_[1] = frame.box_[1][1];
    box_[2] = frame.box_[2][2];

    areaPerLipid_ = box_[0]*box_[1] / numLipids_;
//    printf("APL %8.3f\n", areaPerLipid_);

    makePairs(frame, upperHeads_, lowerHeads_, upperPair_);
    makePairs(frame, lowerHeads_, upperHeads_, lowerPair_);

#pragma omp parallel default(none) shared(frame)
    {
        thicknessWithRef(frame, upperHeads_, lowerHeads_, upperPair_);
        thicknessWithRef(frame, lowerHeads_, upperHeads_, lowerPair_);
    }

    numFrames_++;
}

void Membrane::makePairs(const Frame &frame, const vector<int> &ref,
                         const vector<int> &other, map<int, double> &pairs){
    // For each reference particle in the ref leaflet
    for(const int i : ref){
        double min_dist_2 = box_[2];

        double coords_i[3];
        coords_i[0] = frame.x_[i][0];
        coords_i[1] = frame.x_[i][1];
        coords_i[2] = 0.;
        double coords_j[3];

        // Find the closest reference particle in the other leaflet
        int closest = -1;
        for(const int j : other){
            coords_j[0] = frame.x_[j][0];
            coords_j[1] = frame.x_[j][1];
            coords_j[2] = 0.;

            const double dist_2 = distSqr(coords_i, coords_j);
            if(dist_2 < min_dist_2){
                min_dist_2 = dist_2;
                closest = j;
            }
        }

        coords_i[2] = frame.x_[i][2];
        coords_j[2] = frame.x_[closest][2];
        pairs[i] = fabs(coords_i[2] - coords_j[2]);
    }
}

void Membrane::thicknessWithRef(const Frame &frame, const vector<int> &ref,
                                    const vector<int> &other, const map<int, double> &pairs){
    const double max_box = box_[0] > box_[1] ? box_[0] : box_[1];

    double grid_coords[3];
    double ref_coords[3];

    grid_coords[2] = 0.;
    ref_coords[2] = 0.;

     // For each grid point
    #pragma omp for
    for(int i=0; i<grid_; i++){
        grid_coords[0] = (i + 0.5) * step_[0];

        for(int j=0; j<grid_; j++){
            grid_coords[1] = (j + 0.5) * step_[1];
            double min_dist2 = max_box*max_box;

            // Find closest lipid in reference leaflet
            int closest = -1;
            for(int r : ref){
                ref_coords[0] = frame.x_[r][0];
                ref_coords[1] = frame.x_[r][1];
                const double dist2 = distSqrPlane(grid_coords, ref_coords);
                if(dist2 < min_dist2){
                    closest = r;
                    min_dist2 = dist2;
                }
            }

            // Lookup thickness for the closest lipid
            thickness_(i, j) += pairs.at(closest);
        }
    }
}

double Membrane::mean(){
    return thickness_.mean();
}

void Membrane::normalize(const int smooth_iter){
    // Apply iterations of Gauss-Seidel smoother
    thickness_.smooth(smooth_iter);
    thickness_ /= 2 * numFrames_;
}

void Membrane::printCSV(const std::string &filename, const bool header){
    if(header){
        const string file = filename + ".dat";
        // Backup using small_functions.h
        backup_old_file(file);
        FILE *f = fopen(file.c_str(), "w");

        fprintf(f, "@legend Membrane thickness\n");
        fprintf(f, "@xlabel X (nm)\n");
        fprintf(f, "@ylabel Y (nm)\n");
        fprintf(f, "@xwidth %f\n", box_[0]);
        fprintf(f, "@ywidth %f\n", box_[1]);
        fclose(f);
    }
    thickness_.printCSV(filename);
}

void Membrane::setResolution(const int n){
    grid_ = n;
    thickness_.init(grid_, grid_);
    step_[0] = box_[0] / grid_;
    step_[1] = box_[1] / grid_;
}

void Membrane::reset(){
    thickness_.zero();
    numFrames_ = 0;
}