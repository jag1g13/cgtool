//
// Created by james on 28/04/15.
//

#include "membrane.h"

#include <cmath>

#include "small_functions.h"

using std::string;
using std::vector;

Membrane::Membrane(const string &resname, const string &ref_atom, const int num_atoms, const int num_residues){
    residue_.resname = resname;
    residue_.ref_atom = ref_atom;
    residue_.num_atoms = num_atoms;
    residue_.num_residues = num_residues;
}

//TODO fix this so that ref_atom can be used by name
void Membrane::sortBilayer(const Frame &frame, const int ref_atom){
//    const int ref_atom = frame.nameToNum_.at(residue_.ref_atom);

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

double Membrane::thickness(const Frame &frame){
    vector<double> dists;
    double avg_dist = 0.;

    // For each reference particle in the upper membrane
    for(const int i : upperHeads_){
        double min_dist_2 = frame.box_[2][2];
        int closest;

        double coords_i[3];
        coords_i[0] = frame.atoms_[i].coords[0];
        coords_i[1] = frame.atoms_[i].coords[1];
        coords_i[2] = 0.;
        double coords_j[3];

        // Find the closest reference particle in the lower membrane
        for(const int j : lowerHeads_){
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

        const double dist = sqrt(distSqr(coords_i, coords_j));
        dists.push_back(dist);
        avg_dist += dist;
    }

    avg_dist /= upperHeads_.size();
    return avg_dist;
}
