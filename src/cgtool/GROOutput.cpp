//
// Created by james on 21/08/15.
//

#include "GROOutput.h"

#include <cstdio>
#include <stdexcept>

#include "small_functions.h"

using std::string;

GROOutput::GROOutput(const int natoms, const string &filename){
    natoms_ = natoms;
    if(openFile(filename))
        throw std::runtime_error("ERROR: Could not open GRO file for writing");
}

GROOutput::~GROOutput(){
    closeFile();
}

int GROOutput::openFile(const string &filename){
    backup_old_file(filename);
    file_ = std::fopen(filename.c_str(), "w");
    if(!file_) return 1;

    fprintf(file_, "Generated by CGTOOLS : %s\n", filename.c_str());
    fprintf(file_, "%i\n", natoms_);
    return 0;
}

int GROOutput::closeFile(){
    if(file_) fclose(file_);
    if(!file_) return 0;
    return 1;
}

int GROOutput::writeFrame(const Frame &frame){
    // Print atoms
    Residue &res = frame.residues_[0];
    for(int i=0; i < natoms_; i++){
        const Atom &atom = frame.atoms_[i];
        fprintf(file_, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
                1+(i/res.num_atoms), res.resname.c_str(),
                atom.atom_name.c_str(), i+1,
                atom.coords[0], atom.coords[1], atom.coords[2]);
    }

    // Print box
    fprintf(file_, "%10.5f%10.5f%10.5f\n",
            frame.box_[0][0], frame.box_[1][1], frame.box_[2][2]);

    return 0;
}
