//
// Created by james on 28/04/15.
//

#include "membrane.h"

#include <cmath>

#include "small_functions.h"

using std::string;
using std::vector;

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
    counts_.init(grid_, grid_);
}

void Membrane::sortBilayer(const Frame &frame, const int ref_atom){
//TODO fix this so that ref_atom can be used by name
//    const int ref_atom = frame.nameToNum_.at(residue_.ref_atom);

    // Copy box from Frame - assume orthorhombic
    box_[0] = frame.box_[0][0];
    box_[1] = frame.box_[1][1];
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

double Membrane::thickness(const Frame &frame, const bool with_reset){
    double avg_dist = 0.;

    if(with_reset){
        thickness_.zero();
        counts_.zero();
        numFrames_ = 0;
        avgDist_ = 0.;
    }

    avg_dist += thickness_with_ref(frame, upperHeads_, lowerHeads_);
    avg_dist += thickness_with_ref(frame, lowerHeads_, upperHeads_);
    avg_dist /= 2;
    numFrames_++;
    return avg_dist;
}

double Membrane::thickness_with_ref(const Frame &frame, const vector<int> &ref,
                                    const vector<int> &other){
    double avg_dist = 0.;
    double step[2];
    step[0] = box_[0] / grid_;
    step[1] = box_[1] / grid_;

    // For each reference particle in the ref leaflet
    for(const int i : ref){
        double min_dist_2 = box_[2];
        int closest;

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

        coords_i[2] = frame.atoms_[i].coords[2];
        coords_j[0] = frame.atoms_[closest].coords[0];
        coords_j[1] = frame.atoms_[closest].coords[1];
        coords_j[2] = frame.atoms_[closest].coords[2];

//        const double dist = sqrt(distSqr(coords_i, coords_j));
        const double dist = fabs(coords_i[2] - coords_j[2]);
        avgDist_ += dist;
        avg_dist += dist;

        // Add thickness to grid
        const int x = static_cast<int>((coords_i[0] / step[0]) - 1);
        const int y = static_cast<int>((coords_i[1] / step[1]) - 1);
        thickness_(x, y) += dist;
        counts_(x, y)++;
    }
    return avg_dist / ref.size();
}

double Membrane::mean(){
    thickness_.element_divide(counts_);
    avgDist_ /= (2 * upperHeads_.size() * numFrames_);
    return avgDist_;
}

void Membrane::print_csv(const std::string &filename){
    thickness_.replace_nan();
    thickness_.smooth(2);
    thickness_.print_csv(filename);
}
