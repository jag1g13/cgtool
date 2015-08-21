//
// Created by james on 21/08/15.
//

#ifndef CGTOOL_GROOUTPUT_H
#define CGTOOL_GROOUTPUT_H

#include "trj_output.h"

class GROOutput : public TrjOutput{
protected:
    FILE *file_ = nullptr;
    int openFile(const std::string &filename);
    int closeFile();
public:
    GROOutput(const int natoms, const std::string &filename);
    ~GROOutput();
    int writeFrame(const Frame &frame);

    friend class Frame;
};

#endif //CGTOOL_GROOUTPUT_H
