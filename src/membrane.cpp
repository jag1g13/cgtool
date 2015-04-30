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

Membrane::Membrane(const string &resname, const string &ref_atom,
                   const int num_atoms, const int num_residues,
                   const int grid){
    // Copy arguments into residue_
    //TODO just pass the residue in as a whole
    residue_.resname = resname;
    residue_.ref_atom = ref_atom;
    residue_.num_atoms = num_atoms;
    residue_.num_residues = num_residues;

    // Setup grid
    grid_ = grid;
    thickness_.init(grid_, grid_);
}

void Membrane::sortBilayer(const Frame &frame, const int ref_atom){
//TODO fix this so that ref_atom can be used by name
//    const int ref_atom = frame.nameToNum_.at(residue_.ref_atom);

    // Copy box from Frame - assume orthorhombic
    box_[0] = frame.box_[0][0];
    box_[1] = frame.box_[1][1];    /** Closest in lower leaflet to lipid in upper leaflet */
    std::map<int, int> upperPair_;
    box_[2] = frame.box_[2][2];


    // Calculate average z coord of reference atom
    double avg_z = 0.;
    for(int i=0; i<residue_.num_residues; i++){
        const int num = ref_atom + i * residue_.num_atoms;
        avg_z += frame.atoms_[num].coords[2];
    }
    avg_z /= residue_.num_residues;

    // Separate membrane into upper and lower
    for(int i=0; i<residue_.num_residues; i++){
        const int num = ref_atom + i * residue_.num_atoms;
        const double z = frame.atoms_[num].coords[2];
        if(z < avg_z){
            lowerHeads_.push_back(num);
        }else{
            upperHeads_.push_back(num);
        }
    }
}

void Membrane::thickness(const Frame &frame, const bool with_reset){
    if(with_reset){
        thickness_.zero();
        numFrames_ = 0;
    }

    makePairs(frame, upperHeads_, lowerHeads_, upperPair_);
    makePairs(frame, lowerHeads_, upperHeads_, lowerPair_);

    thicknessWithRef(frame, upperHeads_, lowerHeads_, upperPair_);
    thicknessWithRef(frame, lowerHeads_, upperHeads_, lowerPair_);

    numFrames_++;
}

void Membrane::makePairs(const Frame &frame, const vector<int> &ref,
                         const vector<int> &other, map<int, int> &pairs){
    // For each reference particle in the ref leaflet
    for(const int i : ref){
        double min_dist_2 = box_[2];
        int closest = -1;

        double coords_i[3];
        coords_i[0] = frame.atoms_[i].coords[0];
        coords_i[1] = frame.atoms_[i].coords[1];
        coords_i[2] = 0.;
        double coords_j[3];

        // Find the closest reference particle in the other leaflet
        for(const int j : other){
            coords_j[0] = frame.atoms_[j].coords[0];
            coords_j[1] = frame.atoms_[j].coords[1];
            coords_j[2] = 0.;

            const double dist_2 = distSqr(coords_i, coords_j);
            if(dist_2 < min_dist_2){
                min_dist_2 = dist_2;
                closest = j;
            }
        }

        pairs[i] = closest;
    }
}

void Membrane::thicknessWithRef(const Frame &frame, const vector<int> &ref,
                                    const vector<int> &other, const map<int, int> &pairs){
    double step[2];
    step[0] = box_[0] / grid_;
    step[1] = box_[1] / grid_;
    const double max_box = box_[0] > box_[1] ? box_[0] : box_[1];

    double grid_coords[3];
    double ref_coords[3];
    double other_coords[3];

    grid_coords[2] = 0.;
    ref_coords[2] = 0.;

    // For each grid point
    for(int i=0; i<grid_; i++){
        grid_coords[0] = i * step[0];

        for(int j=0; j<grid_; j++){
            grid_coords[1] = j * step[1];
            double min_dist2 = max_box*max_box;

            // Find closest in reference leaflet
            int closest = -1;
            for(int r : ref){
                ref_coords[0] = frame.atoms_[r].coords[0];
                ref_coords[1] = frame.atoms_[r].coords[1];
                const double dist2 = distSqr(grid_coords, ref_coords);
                if(dist2 < min_dist2){
                    closest = r;
                    min_dist2 = dist2;
                }
            }

            const int other_r = pairs.at(closest);
            ref_coords[2] = frame.atoms_[closest].coords[2];
            other_coords[2] = frame.atoms_[other_r].coords[2];
            thickness_(i, j) += fabs(ref_coords[2] - other_coords[2]);
        }
    }
}

double Membrane::mean(){
    return thickness_.mean();
}

void Membrane::normalize(){
    thickness_ /= 2 * numFrames_;
}

void Membrane::printCSV(const std::string &filename){
    thickness_.smooth(2);
    thickness_.print_csv(filename);
}
