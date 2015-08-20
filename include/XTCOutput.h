//
// Created by james on 20/08/15.
//

#ifndef CGTOOL_XTCOUTPUT_H
#define CGTOOL_XTCOUTPUT_H

#include "trj_output.h"

class XTCOutput : public TrjOutput{
protected:
    rvec *x_ = nullptr;
    float box_[3][3];
    XDRFILE *file_;
    int openFile(const std::string &filename);
    int closeFile();
public:
    XTCOutput(const int natoms);
    ~XTCOutput();
    int writeFrame(const Frame &frame);

    friend class Frame;
};

#endif //CGTOOL_XTCOUTPUT_H
