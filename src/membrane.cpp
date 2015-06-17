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
    residuePPL_.resize(residues_.size());

    prepCSVAreaPerLipid();

    curvMean_.alloc(grid_, grid_);
    curvGaussian_.alloc(grid_, grid_);
}

void Membrane::sortBilayer(const Frame &frame, const int blocks){
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
        printf("Residue %6s ref %6s%5i\n", res.resname.c_str(), res.ref_atom_name.c_str(), res.ref_atom);
        numLipids_ += res.num_residues;
        for(int i = 0; i < res.num_residues; i++){
            const int num = res.ref_atom + i * res.num_atoms + res.start;
            const int x = std::max(int(frame.x_[num][0] * blocks / box_[0]), 0);
            const int y = std::min(int(frame.x_[num][1] * blocks / box_[1]), blocks-1);
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
            const int x = std::max(int(frame.x_[num][0] * blocks / box_[0]), 0);
            const int y = std::min(int(frame.x_[num][1] * blocks / box_[1]), blocks-1);
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

    makePairs(frame, upperHeads_, lowerHeads_, upperPair_);
    makePairs(frame, lowerHeads_, upperHeads_, lowerPair_);

#pragma omp parallel default(none) shared(frame)
    {
        closestLipid(frame, upperHeads_, upperPair_, closestUpper_);
        closestLipid(frame, lowerHeads_, lowerPair_, closestLower_);
    }

    areaPerLipid(closestUpper_);
    areaPerLipid(closestLower_);
    printCSVAreaPerLipid(frame.time_);

    curvature(closestUpper_, closestLower_, frame);
//    printCSVCurvature("curvature_" + std::to_string(frame.num_));
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

            const double dist_2 = distSqrPlane(coords_i, coords_j);
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

void Membrane::closestLipid(const Frame &frame, const std::vector<int> &ref,
                            const std::map<int, double> &pairs, LightArray<int> &closest){
    const double max_box = box_[0] > box_[1] ? box_[0] : box_[1];

    double grid_coords[3];
    double ref_coords[3];

    grid_coords[2] = 0.;
    ref_coords[2] = 0.;

#pragma omp for
    for(int i=0; i<grid_; i++){
        grid_coords[0] = (i + 0.5) * step_[0];

        for(int j=0; j<grid_; j++){
            grid_coords[1] = (j + 0.5) * step_[1];
            double min_dist2 = max_box*max_box;

            // Find closest lipid in reference leaflet
            closest(i, j) = -1;
            for(int r : ref){
                ref_coords[0] = frame.x_[r][0];
                ref_coords[1] = frame.x_[r][1];
                const double dist2 = distSqrPlane(grid_coords, ref_coords);
                if(dist2 < min_dist2){
                    closest(i, j) = r;
                    min_dist2 = dist2;
                }
            }

            thickness_(i, j) += pairs.at(closest.at(i, j));
        }
    }
}

void Membrane::areaPerLipid(const LightArray<int> &closest){
    for(int &n : residuePPL_) n = 0;

    for(int i=0; i<grid_; i++){
        for(int j=0; j<grid_; j++){
            const int lipid = closest.at(i, j);
            for(int k=0; k<residues_.size(); k++){
                if(lipid >= residues_[k].start && lipid < residues_[k].end){
                    residuePPL_[k]++;
                }
            }
        }
    }

}

void Membrane::curvature(const LightArray<int> &upper, const LightArray<int> &lower,
                         const Frame &frame){
    LightArray<double> avg_z(grid_, grid_);

    // Calculate average z coord on grid
    for(int i=0; i<grid_; i++){
        for(int j=0; j<grid_; j++){
            avg_z(i, j) = (frame.x_[upper.at(i, j)][2] + frame.x_[lower.at(i, j)][2]) / 2.;
        }
    }

    LightArray<double> respect_to_x(grid_, grid_);
    LightArray<double> respect_to_y(grid_, grid_);

    const double inv_h2_x = 1. / (step_[0]*step_[0]);
    const double inv_h2_y = 1. / (step_[1]*step_[1]);

    double curv_x_avg = 0.;
    double curv_y_avg = 0.;

    // Do finite differences wrt x and y
    for(int i=1; i<grid_-1; i++){
        for(int j=1; j<grid_-1; j++){
            respect_to_x(i, j) = inv_h2_x * (avg_z(i+1, j) + avg_z(i-1, j) - 2*avg_z(i, j));
            respect_to_y(i, j) = inv_h2_y * (avg_z(i, j+1) + avg_z(i, j-1) - 2*avg_z(i, j));

            curv_x_avg += respect_to_x(i, j);
            curv_y_avg += respect_to_y(i, j);
        }
    }

    curv_x_avg /= grid_*grid_;
    curv_y_avg /= grid_*grid_;

//    printf("Calculate\n");
    for(int i=0; i<grid_; i++){
        for(int j=0; j<grid_; j++){
            curvMean_(i, j) = (respect_to_x.at(i, j) + respect_to_y.at(i, j)) / 2.;
            curvGaussian_(i, j) = respect_to_x.at(i, j) * respect_to_y.at(i, j);
//            printf("%8.3f", curvMean_.at(i, j));
        }
//        printf("\n");
    }

//    printf("%8.3f%8.3f\n", curv_x_avg, curv_y_avg);
}

void Membrane::printCSVCurvature(const std::string &filename) const{
    const string file = filename + ".dat";
    // Backup using small_functions.h
    backup_old_file(file);
    FILE *f = fopen(file.c_str(), "w");
    if(f == nullptr) throw std::runtime_error("Could not open output file.");

    if(header_){
        fprintf(f, "@legend Membrane mean curvature\n");
        fprintf(f, "@xlabel X (nm)\n");
        fprintf(f, "@ylabel Y (nm)\n");
        fprintf(f, "@xwidth %f\n", box_[0]);
        fprintf(f, "@ywidth %f\n", box_[1]);
    }

//    printf("Printing\n");
    for(int i=0; i<grid_; i++){
        for(int j=0; j<grid_; j++){
            fprintf(f, "%8.3f", curvMean_.at(i, j));
//            printf("%8.3f", curvMean_.at(i, j));
        }
        fprintf(f, "\n");
//        printf("\n");
    }

    fclose(f);
}

void Membrane::prepCSVAreaPerLipid(){
    const string file = "APL.dat";
    // Backup using small_functions.h
    backup_old_file(file);
    aplFile_ = fopen(file.c_str(), "w");

    if(header_){
        fprintf(aplFile_, "@legend Area Per Lipid\n");
        fprintf(aplFile_, "@xlabel Time (ps)\n");
        fprintf(aplFile_, "@ylabel APL (nm^2)\n");
    }
}

void Membrane::printCSVAreaPerLipid(const float time) const{
    fprintf(aplFile_, "%12.3f", time);
    for(int i=0; i<residuePPL_.size(); i++){
        const int num = residuePPL_[i];
        const double area = num * step_[0] * step_[1];
        const double APL = area / residues_[i].num_residues;
        fprintf(aplFile_, "%8.3f", APL);
    }
    fprintf(aplFile_, "\n");
}

double Membrane::mean() const{
    return thickness_.mean();
}

void Membrane::normalize(const int smooth_iter){
    // Apply iterations of Gauss-Seidel smoother
    thickness_.smooth(smooth_iter);
    thickness_ /= 2 * numFrames_;
}

void Membrane::printCSV(const std::string &filename) const{
    if(header_){
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
    // Print CSV - true suppresses backup - file has been opened already
    thickness_.printCSV(filename, true);
}

void Membrane::setResolution(const int n){
    grid_ = n;
    thickness_.init(grid_, grid_);
    step_[0] = box_[0] / grid_;
    step_[1] = box_[1] / grid_;

    closestUpper_.alloc(n, n);
    closestLower_.alloc(n, n);
}

void Membrane::reset(){
    thickness_.zero();
    numFrames_ = 0;
}