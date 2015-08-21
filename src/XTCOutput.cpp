//
// Created by james on 20/08/15.
//

#include "XTCOutput.h"

#include "xdrfile_xtc.h"

#include "small_functions.h"

using std::string;

XTCOutput::XTCOutput(const int natoms, const string &filename){
    natoms_ = natoms;
    x_ = new rvec[natoms_];
    openFile(filename);
}

XTCOutput::~XTCOutput(){
    closeFile();
    delete[] x_;
}

int XTCOutput::openFile(const string &filename){
    char mode[2] = {'r', 'w'};
    backup_old_file(filename);
    file_ = xdrfile_open(filename.c_str(), &mode[1]);
    if(file_) return 0;
    return 1;
}

int XTCOutput::closeFile(){
    if(file_) xdrfile_close(file_);
    if(!file_) return 0;
    return 1;
}

int XTCOutput::writeFrame(const Frame &frame){
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