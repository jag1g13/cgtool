//
// Created by james on 21/08/15.
//

#include "GROOutput.h"
#include "small_functions.h"

using std::string;

GROOutput::GROOutput(const int natoms, const string &filename){
    natoms_ = natoms;
    if(openFile(filename))
        throw std::runtime_error("ERROR: Could not open GRO file for writing");
    std::stringstream stream;
    stream << "Generated by CGTOOL : " << (*residues_)[0].resname << "\n" << natoms_ << "\n";
    fprintf(gro, "%s", stream.str().c_str());

    fprintf("")
}

GROOutput::~GROOutput(){
    closeFile();
    delete[] x_;
}

int GROOutput::openFile(const string &filename){
    backup_old_file(filename);
    file_ = std::fopen(filename.c_str(), "w");
    if(file_) return 0;
    return 1;
}

int GROOutput::closeFile(){
    if(file_) xdrfile_close(file_);
    if(!file_) return 0;
    return 1;
}

int GROOutput::writeFrame(const Frame &frame){
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            box_[i][j] = frame.box_[i][j];
        }
    }

    for(int i=0; i<natoms_; i++){
        x_[i][0] = float(frame.atoms_[i].coords[0]);
        x_[i][1] = float(frame.atoms_[i].coords[1]);
        x_[i][2] = float(frame.atoms_[i].coords[2]);
    }

    return exdrOK == write_xtc(file_, natoms_, frame.step_,
                               frame.time_, box_, x_, frame.prec_);
}
