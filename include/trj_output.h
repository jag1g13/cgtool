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
    /** \brief Number of atoms to be written to output. */
    int natoms_;

    /** \brief Open and prepare output file. */
    virtual int openFile(const std::string &filename) = 0;
    /** \brief Close output file. */
    virtual int closeFile() = 0;

public:
    /** \brief Write a Frame to output file.  Pure virtual function. */
    virtual int writeFrame(const Frame &frame) = 0;

    /** \brief Empty destructor to be overwritten. */
    virtual ~TrjOutput(){};

    friend class Frame;
};


#endif //CGTOOL_TRJOUTPUT_H
