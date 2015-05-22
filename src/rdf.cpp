//
// Created by james on 22/05/15.
//

#include "rdf.h"

using std::vector;

RDF::RDF(const vector<Residue> &residues){
    residues_ = residues;
}

void RDF::calculateRDF(const Frame &frame){
    for(const Residue &res : residues_){
        for(int i=0; i<res.num_residues; i++){
            // Get number of ref atom for this residue
            const int atom = res.start + i*res.num_atoms + res.ref_atom;

            double ref_coords[3];
            ref_coords[0] = frame.x_[atom][0];
            ref_coords[1] = frame.x_[atom][1];
            ref_coords[2] = frame.x_[atom][2];
        }
    }
}