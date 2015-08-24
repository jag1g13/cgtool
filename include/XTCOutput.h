//
// Created by james on 20/08/15.
//

#ifndef CGTOOL_XTCOUTPUT_H
#define CGTOOL_XTCOUTPUT_H

#include "trj_output.h"
#include "xdrfile.h"

class XTCOutput : public TrjOutput{
protected:
    /** \brief Stores atoms copied from Frame. */
    rvec *x_ = nullptr;
    /** \brief Stores box vectors copied from Frame. */
    float box_[3][3];
    /** \brief Output file handle. */
    XDRFILE *file_;

    /** \brief Open and prepare output file. */
    int openFile(const std::string &filename);
    /** \brief Close output file. */
    int closeFile();
public:
    /** \brief Constructor.  Calls openFile(). */
    XTCOutput(const int natoms, const std::string &filename);
    /** \brief Destructor.  Calls closeFile(). */
    ~XTCOutput();

    /** \brief Write a Frame to output file. */
    int writeFrame(const Frame &frame);

    friend class Frame;
};

#endif //CGTOOL_XTCOUTPUT_H
