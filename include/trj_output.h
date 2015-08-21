//
// Created by james on 20/08/15.
//

#ifndef CGTOOL_TRJOUTPUT_H
#define CGTOOL_TRJOUTPUT_H

#include <string>
#include <stdexcept>

#include "frame.h"
//class Frame;

class TrjOutput{
protected:
    int natoms_;
    virtual int openFile(const std::string &filename) = 0;
    virtual int closeFile() = 0;

public:
    virtual int writeFrame(const Frame &frame) = 0;
    virtual ~TrjOutput(){};

    friend class Frame;
};


#endif //CGTOOL_TRJOUTPUT_H
