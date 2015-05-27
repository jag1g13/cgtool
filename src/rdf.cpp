//
// Created by james on 22/05/15.
//

#include "rdf.h"

#include <cmath>

#include "small_functions.h"

using std::vector;

RDF::RDF(const vector<Residue> &residues){
    init(residues, 2., 100);
}

void RDF::init(const vector<Residue> &residues, const double cutoff, const int resolution){
    residues_ = residues;
    cutoff_ = cutoff;
    resolution_ = resolution;
    grid_ = static_cast<int>(cutoff * resolution_);
    histogram_.init(grid_);
    rdf_.init(grid_);
}

void RDF::calculateRDF(const Frame &frame){
    // Calculate average number density in cell
    const double volume = frame.box_[0][0] * frame.box_[1][1] * frame.box_[2][2];
    const double density = frame.residues_[0].num_residues / volume;

    Residue &res = residues_[0];
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

            const double dist2 = distSqr(R_a, R_b);
            // Set cutoff at 2nm
            if(dist2 > cutoff_*cutoff_) continue;
            const double dist = sqrt(dist2);
            const int loc = static_cast<int>(dist * resolution_);
            histogram_.increment(loc);
        }
    }

    // Populate rdf_ with reciprocal of expected number per shell
    const double prefactor = (4. / 3.) * M_PI;
    for(int i=0; i<grid_; i++){
        const double r_inner = i * resolution_;
        const double r_outer = r_inner + resolution_;
        const double v_inner = r_inner * r_inner * r_inner;
        const double v_outer = r_outer * r_outer * r_outer;
        rdf_(i) = density * prefactor * (v_outer - v_inner);
    }

    rdf_.elementMultiply(histogram_);

    histogram_.print();
    rdf_.print();
//    rdf_.printCSV("rdf");
}