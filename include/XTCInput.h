//
// Created by james on 24/08/15.
//

#ifndef CGTOOL_XTCINPUT_H
#define CGTOOL_XTCINPUT_H

#include "TrjInput.h"
#include "xdrfile.h"

class XTCInput : TrjInput{
protected:
    /** \brief Stores atoms from XTC. */
    rvec *x_ = nullptr;
    /** \brief Stores box vectors from XTC. */
    float box_[3][3];
    /** \brief Input file handle. */
    XDRFILE *file_ = nullptr;
    /** \brief Current timestep */
    float step_;
    /** \brief Current time */
    float time_;
    /** \brief XTC precision */
    float prec_;

    /** \brief Open and prepare input file. */
    int openFile(const std::string &filename);
    /** \brief Close input file. */
    int closeFile();
public:
    /** \brief Constructor.  Calls openFile(). */
    XTCInput(const int natoms, const std::string &filename);
    /** \brief Destructor.  Calls closeFile(). */
    ~XTCInput();

    /** \brief Read a Frame from input file. */
    int readFrame(Frame &frame);

    friend class Frame;
};


#endif //CGTOOL_XTCINPUT_H
