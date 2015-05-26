//
// Created by james on 22/05/15.
//

#include "rdf.h"

#include <cmath>

#include "small_functions.h"

using std::vector;

RDF::RDF(const vector<Residue> &residues){
    residues_ = residues;
}

void RDF::calculateRDF(const Frame &frame){
    for(const Residue &res : residues_){
        for(int i=0; i<res.num_residues; i++){
            // Get number of ref atom for this residue
            const int atom_a = res.start + i*res.num_atoms + res.ref_atom;

            double R_a[3];
            R_a[0] = frame.x_[atom_a][0];
            R_a[1] = frame.x_[atom_a][1];
            R_a[2] = frame.x_[atom_a][2];

            for(int j=0; j<res.num_residues; j++){
                if(i == j) continue;

                const int atom_b = res.start + i*res.num_atoms + res.ref_atom;
                double R_b[3];
                R_b[0] = frame.x_[atom_b][0];
                R_b[1] = frame.x_[atom_b][1];
                R_b[2] = frame.x_[atom_b][2];

                const double dist = sqrt(distSqr(R_a, R_b));
                // Convert distances to 0.1A - usual resolution of XTC files
                const int loc = int(dist * 100);
                histogram_.increment(loc);
            }
        }
    }
}