//
// Created by james on 24/08/15.
//

#include "XTCInput.h"

#include <cstdio>

#include "xdrfile_xtc.h"

using std::string;
using std::printf;

XTCInput::XTCInput(const string &filename){
    // How many atoms?  Prepare Frame for reading
    int status = read_xtc_natoms(filename.c_str(), &natoms_);
    if(status != exdrOK) throw std::runtime_error("Could not open input XTC for reading");

    x_ = new rvec[natoms_];
    if(openFile(filename)) throw std::runtime_error("Error reading initial frame from XTC");
}

XTCInput::~XTCInput(){
    closeFile();
    if(x_) delete[] x_;
}

int XTCInput::openFile(const std::string &filename){
    file_ = xdrfile_open(filename.c_str(), "r");
    int status = read_xtc(file_, natoms_, &step_, &time_, box_, x_, &prec_);
    if(status != exdrOK) return 1;

    // Check box vectors
    bool cubic = true;
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            if(i!=j && box_[i][j] > 0) cubic = false;
        }
    }
    if(!cubic) printf("NOTE: Input box is not cubic\n");

    return 0;
}

int XTCInput::closeFile(){
    if(file_) xdrfile_close(file_);
    return 0;
}

int XTCInput::readFrame(Frame &frame){
    int status = read_xtc(file_, natoms_, &step_, &time_, box_, x_, &prec_);
    if(status != exdrOK) return 1;

//    pbcAtom();
//    copyCoordsIntoAtoms();

    // Copy data into Frame
    frame.step_ = step_;
    frame.time_ = time_;
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            frame.box_[i][j] = box_[i][j];
        }
    }
    for(int i=0; i<frame.numAtoms_ && i<natoms_; i++){
        frame.atoms_[i].coords[0] = x_[i][0];
        frame.atoms_[i].coords[1] = x_[i][1];
        frame.atoms_[i].coords[2] = x_[i][2];
    }

    return 0;
}