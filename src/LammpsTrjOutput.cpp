//
// Created by james on 27/08/15.
//

#include "LammpsTrjOutput.h"

#include <cstdio>

#include "small_functions.h"

using std::string;

LammpsTrjOutput::LammpsTrjOutput(const int natoms, const string &filename){
    natoms_ = natoms;
    if(openFile(filename))
        throw std::runtime_error("ERROR: Could not open Lammps Trj file for writing");
}

LammpsTrjOutput::~LammpsTrjOutput(){
    closeFile();
}

int LammpsTrjOutput::openFile(const string &filename){
    backup_old_file(filename);
    file_ = std::fopen(filename.c_str(), "w");
    if(!file_) return 1;
    return 0;
}

int LammpsTrjOutput::closeFile(){
    if(file_) fclose(file_);
    if(!file_) return 0;
    return 1;
}

int LammpsTrjOutput::writeFrame(const Frame &frame){
    // Have to multiply all coords by 10 - Gromacs is in A, Lammps in nm
    double box[3];
    for(int i=0; i<3; i++){
        box[i] = 10 * frame.box_[i][i] / 2;
    }

    // Print headers
    fprintf(file_, "ITEM: TIMESTEP\n%d\n", frame.num_);
    fprintf(file_, "ITEM: NUMBER OF ATOMS\n%d\n", frame.numAtoms_);
    fprintf(file_, "ITEM: BOX BOUNDS pp pp pp\n%f %f\n%f %f\n%f %f\n",
            -box[0], box[0], -box[1], box[1], -box[2], box[2]);
    fprintf(file_, "ITEM: ATOMS id type mol x y z mux muy muz mass diameter\n");

    for(int i=0; i < natoms_; i++){
        const Atom &atom = frame.atoms_[i];
        fprintf(file_, " %6d %4d %4d %10.4f %10.4f %10.4f %9.5f %9.5f %9.5f %4.1f %4.1f\n",
                i+1, i+1, 1,
                10*atom.coords[0]-box[0], 10*atom.coords[1]-box[1], 10*atom.coords[2]-box[2],
                10*atom.dipole[0], 10*atom.dipole[1], 10*atom.dipole[2],
                atom.mass, 2.7f);
    }

    return 0;
}
