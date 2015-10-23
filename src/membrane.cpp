//
// Created by james on 28/04/15.
//

#include "membrane.h"

#include <cmath>
#include <cstdio>
#include <array>
#include <boost/algorithm/clamp.hpp>

#include "small_functions.h"

using std::string;
using std::vector;
using std::array;
using std::map;
using std::set;
using boost::algorithm::clamp;

Membrane::Membrane(const vector<Residue> &residues, const Frame &frame,
                   const int resolution, const int blocks, const bool header) :
        residues_(residues), header_(header){
    setResolution(resolution);
    sortBilayer(frame, blocks);
    prepCSVAreaPerLipid();
    prepCSVAvgThickness();
}

Membrane::~Membrane(){
    fclose(aplFile_);
    fclose(avgFile_);
};

void Membrane::sortBilayer(const Frame &frame, const int blocks){
    // Copy box from Frame - assume orthorhombic
    box_[0] = frame.box_[0][0];
    box_[1] = frame.box_[1][1];
    box_[2] = frame.box_[2][2];

    // Find middle of membrane in z coord
    // Do this in blocks to account for curvature - more curvature needs more blocks
    LightArray<double> block_avg_z(blocks, blocks);
    LightArray<int> block_tot_residues(blocks, blocks);

    for(const Residue &res : residues_){
        if(res.ref_atom < 0) continue;
        numLipids_ += res.num_residues;
        for(int i = 0; i < res.num_residues; i++){
            const int num = res.ref_atom + i * res.num_atoms + res.start;
            const int x = int(frame.atoms_[num].coords[0] * blocks / box_[0]) % blocks;
            const int y = int(frame.atoms_[num].coords[1] * blocks / box_[1]) % blocks;
            block_avg_z(x, y) += frame.atoms_[num].coords[2];
            block_tot_residues(x, y)++;
        }
    }

    for(int i=0; i<blocks; i++){
        for(int j=0; j<blocks; j++){
            block_avg_z(i, j) /= block_tot_residues(i, j);
        }
    }

    // Separate membrane into upper and lower
    for(const Residue &res : residues_){
        if(res.ref_atom < 0) continue;
        int num_in_leaflet[2] = {0, 0};
        for(int i = 0; i < res.num_residues; i++){
            const int num = res.ref_atom + i * res.num_atoms + res.start;
            const int x = int(frame.atoms_[num].coords[0] * blocks / box_[0]) % blocks;
            const int y = int(frame.atoms_[num].coords[1] * blocks / box_[1]) % blocks;
            const double z = frame.atoms_[num].coords[2];

            if(z < block_avg_z(x, y)){
                lowerHeads_.insert(num);
                num_in_leaflet[0]++;
                lowerNumRes_[res.resname]++;
            }else{
                upperHeads_.insert(num);
                num_in_leaflet[1]++;
                upperNumRes_[res.resname]++;
            }
        }
        printf("%5s: %'4d lower, %'4d upper\n",
               res.resname.c_str(), num_in_leaflet[0], num_in_leaflet[1]);
    }
}

double Membrane::thickness(const Frame &frame, const bool with_reset){
    if(with_reset) reset();

    // Copy box from Frame - assume orthorhombic
    box_[0] = frame.box_[0][0];
    box_[1] = frame.box_[1][1];
    box_[2] = frame.box_[2][2];
    step_[0] = box_[0] / grid_;
    step_[1] = box_[1] / grid_;

    makePairs(frame, upperHeads_, lowerHeads_, upperPair_);
    makePairs(frame, lowerHeads_, upperHeads_, lowerPair_);

    double avg_thickness = 0.;
    avg_thickness += closestLipid(frame, upperHeads_, upperPair_, upperResPPL_, closestUpper_);
    avg_thickness += closestLipid(frame, lowerHeads_, lowerPair_, lowerResPPL_, closestLower_);
    avg_thickness /= 2 * grid_ * grid_;
    fprintf(avgFile_, "%8.3f%8.3f\n", frame.time_, avg_thickness);

    numFrames_++;
    return avg_thickness;
}

void Membrane::makePairs(const Frame &frame, const set<int> &ref,
                         const set<int> &other, map<int, double> &pairs){
    // For each reference particle in the ref leaflet
    for(const int i : ref){
        double min_dist_2 = box_[0] * box_[1];

        array<double, 3> coords_i, coords_j;
        coords_i[0] = frame.atoms_[i].coords[0];
        coords_i[1] = frame.atoms_[i].coords[1];

        // Find the closest reference particle in the other leaflet
        int closest = -1;
        for(const int j : other){
            coords_j[0] = frame.atoms_[j].coords[0];
            coords_j[1] = frame.atoms_[j].coords[1];

            const double dist_2 = distSqrPlane(coords_i, coords_j);
            if(dist_2 < min_dist_2){
                min_dist_2 = dist_2;
                closest = j;
            }
        }

        if(closest == -1) throw std::logic_error("Could not find closest partner lipid in membrane");

        coords_i[2] = frame.atoms_[i].coords[2];
        coords_j[2] = frame.atoms_[closest].coords[2];
        pairs[i] = fabs(coords_i[2] - coords_j[2]);
    }
}

//TODO optimise this - it takes >80% of runtime
double Membrane::closestLipid(const Frame &frame, const set<int> &ref,
                              const map<int, double> &pairs,
                              map<string, int> &resPPL, LightArray<int> &closest){
    const double max_box = box_[0] > box_[1] ? box_[0] : box_[1];

    array<double, 2> grid_coords;
    double sum = 0.;

    const size_t ref_len = ref.size();
    vector<array<double, 2>> ref_cache(ref_len);
    vector<int> ref_lookup(ref_len);
    int it = 0;
    for(const int r : ref){
        ref_cache[it][0] = frame.atoms_[r].coords[0];
        ref_cache[it][1] = frame.atoms_[r].coords[1];
        ref_lookup[it] = r;
        it++;
    }

#pragma omp parallel for default(none) \
 shared(frame, ref, pairs, closest, ref_cache, ref_lookup, resPPL) \
 private(grid_coords) \
 reduction(+: sum)
    for(int i=0; i<grid_; i++){
        grid_coords[0] = (i + 0.5) * step_[0];

        for(int j=0; j<grid_; j++){
            grid_coords[1] = (j + 0.5) * step_[1];
            double min_dist2 = 2*max_box*max_box;

            // Find closest lipid in reference leaflet
            int closest_int = -1;
            for(int k=0; k<ref_len; k++){
                const double dist2 = distSqrPlane(grid_coords, ref_cache[k]);
                closest_int = dist2 < min_dist2 ? k : closest_int;
                min_dist2 = std::min(min_dist2, dist2);
            }

            const int close_ref = ref_lookup[closest_int];
            closest(i, j) = close_ref;
            const double tmp = pairs.at(close_ref);
            for(const Residue &res : residues_){
                if(close_ref >= res.start && close_ref < res.end){
#pragma omp atomic update
                    resPPL[res.resname]++;
                }
            }

            sum += tmp;
            thickness_(i, j) += tmp;
        }
    }

    return sum;
}


void Membrane::curvature(const Frame &frame){
    LightArray<double> avg_z(grid_, grid_);

    LightArray<double> respect_to_x(grid_, grid_);
    LightArray<double> respect_to_y(grid_, grid_);

    const double inv_h2_x = 1. / (step_[0] * step_[0]);
    const double inv_h2_y = 1. / (step_[1] * step_[1]);

    double curv_x_avg = 0.;
    double curv_y_avg = 0.;

    // Calculate average z coord on grid
    for (int i = 0; i < grid_; i++) {
        for (int j = 0; j < grid_; j++) {
            avg_z(i, j) = (frame.atoms_[closestUpper_.at(i, j)].coords[2] +
                           frame.atoms_[closestLower_.at(i, j)].coords[2]) / 2.;
        }
    }

    // Do finite differences wrt x and y
    for (int i = 1; i < grid_ - 1; i++) {
        for (int j = 1; j < grid_ - 1; j++) {
            respect_to_x(i, j) = inv_h2_x * (avg_z(i + 1, j) + avg_z(i - 1, j) - 2 * avg_z(i, j));
            respect_to_y(i, j) = inv_h2_y * (avg_z(i, j + 1) + avg_z(i, j - 1) - 2 * avg_z(i, j));

            curv_x_avg += respect_to_x(i, j);
            curv_y_avg += respect_to_y(i, j);
        }
    }

    curv_x_avg /= grid_*grid_;
    curv_y_avg /= grid_*grid_;

    for(int i=0; i<grid_; i++){
        for(int j=0; j<grid_; j++){
            curvMean_(i, j) = (respect_to_x.at(i, j) + respect_to_y.at(i, j)) / 2.;
            curvGaussian_(i, j) = respect_to_x.at(i, j) * respect_to_y.at(i, j);
        }
    }

    curvMean_.smooth(5);
    curvGaussian_.smooth(5);
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
    if(!aplFile_) throw std::runtime_error("Could not open output file");

    if(header_){
        fprintf(aplFile_, "@legend Area Per Lipid\n");
        fprintf(aplFile_, "@xlabel Time (ps)\n");
        fprintf(aplFile_, "@ylabel APL (nm^2)\n");
    }

    for(const auto &ppl : upperResPPL_) fprintf(aplFile_, "%5s", ppl.first.c_str());
    for(const auto &ppl : lowerResPPL_) fprintf(aplFile_, "%5s", ppl.first.c_str());
    fprintf(aplFile_, "\n");
}

void Membrane::printCSVAreaPerLipid(const float time) const{
    fprintf(aplFile_, "%12.3f", time);
    for(const auto &ppl : upperResPPL_){
        const int num = ppl.second;
        const double area = num * step_[0] * step_[1];
        const double APL = area / (numFrames_ * upperNumRes_.at(ppl.first));
        fprintf(aplFile_, "%12.3f", APL);
    }
    for(const auto &ppl : lowerResPPL_){
        const int num = ppl.second;
        const double area = num * step_[0] * step_[1];
        const double APL = area / (numFrames_ * lowerNumRes_.at(ppl.first));
        fprintf(aplFile_, "%12.3f", APL);
    }
    fprintf(aplFile_, "\n");
}

void Membrane::prepCSVAvgThickness(){
    const string file = "avg_thickness.dat";
    // Backup using small_functions.h
    backup_old_file(file);
    avgFile_ = fopen(file.c_str(), "w");
    if(!avgFile_) throw std::runtime_error("Could not open output file");

    if(header_){
        fprintf(aplFile_, "@legend Average Thickness\n");
        fprintf(aplFile_, "@xlabel Time (ps)\n");
        fprintf(aplFile_, "@ylabel Thickness (nm)\n");
    }
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

        // Print CSV - true suppresses backup - file has been opened already
        thickness_.printCSV(filename, true);
    }else{
        thickness_.printCSV(filename);
    }
}

void Membrane::setResolution(const int n){
    grid_ = n;
    thickness_.alloc(grid_, grid_);

    closestUpper_.alloc(n, n);
    closestLower_.alloc(n, n);
    curvMean_.alloc(grid_, grid_);
    curvGaussian_.alloc(grid_, grid_);
}

void Membrane::reset(){
    thickness_.zero();
    for(auto &ppl : upperResPPL_) ppl.second = 0;
    for(auto &ppl : lowerResPPL_) ppl.second = 0;
    numFrames_ = 0;
}
