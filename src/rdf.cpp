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
    const double box[3] = {frame.box_[0][0], frame.box_[1][1], frame.box_[2][2]};
    const double volume = box[0] * box[1] * box[2];
    density_ += frame.residues_[0].num_residues / volume;

    Residue &res = residues_[0];
    #pragma omp parallel for default(none) shared(res, frame, box)
    for(int i=0; i<res.num_residues; i++){
        // Get number of ref atom for this residue
        const int atom_a = res.start + i*res.num_atoms + res.ref_atom;

        double R_a[3];
        R_a[0] = frame.x_[atom_a][0];
        R_a[1] = frame.x_[atom_a][1];
        R_a[2] = frame.x_[atom_a][2];

        int pbc_axis[3];

        // Assumes cubic/orthorhombic box
        for(int ii=0; ii<3; ii++){
            pbc_axis[ii] = 0;
            if(R_a[ii] < cutoff_) pbc_axis[ii] = -1;
            if(R_a[ii] > box[ii] - cutoff_) pbc_axis[ii] = 1;
        }

        for(int j=0; j<res.num_residues; j++){
            if(i == j) continue;

            const int atom_b = res.start + j*res.num_atoms + res.ref_atom;
            double R_b[3];
            R_b[0] = frame.x_[atom_b][0];
            R_b[1] = frame.x_[atom_b][1];
            R_b[2] = frame.x_[atom_b][2];

            // Loop over centre and neighbouring boxes
            for(int ii=-1; ii<=1; ii++){
                if(pbc_axis[0] != ii && ii != 0) continue;
                for(int jj=-1; jj<=1; jj++){
                    if(pbc_axis[1] != jj && jj != 0) continue;
                    for(int kk=-1; kk<=1; kk++){
                        if(pbc_axis[2] != kk && kk != 0) continue;

                        double R_b_adj[3];
                        R_b_adj[0] = R_b[0] + ii*box[0];
                        R_b_adj[1] = R_b[1] + jj*box[1];
                        R_b_adj[2] = R_b[2] + kk*box[2];

                        const double dist2 = distSqr(R_a, R_b_adj);
                        if(dist2 > cutoff_*cutoff_) continue;
                        const double dist = sqrt(dist2);
                        const int loc = static_cast<int>(dist * resolution_);
                        histogram_.increment(loc);
                    }
                }
            }
        }
    }

    frames_++;
}

void RDF::normalize(){
    const double recip_frames = 1. / frames_;
    density_ *= recip_frames;

    // Populate rdf_ with reciprocal of expected number per shell
    const double prefactor = (4. / 3.) * M_PI;
    const double r_scale = cutoff_ / resolution_;
    for(int i=0; i<grid_; i++){
        const double r_inner = i * r_scale;
        const double r_outer = r_inner + r_scale;
        const double v_inner = r_inner * r_inner * r_inner;
        const double v_outer = r_outer * r_outer * r_outer;
        rdf_(i) = 1. / (density_ * prefactor * (v_outer - v_inner));
    }

    rdf_.elementMultiply(histogram_);
    rdf_ *= recip_frames;
    rdf_.printCSV("rdf");
    frames_ = 0;
}