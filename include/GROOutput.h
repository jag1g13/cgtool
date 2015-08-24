//
// Created by james on 21/08/15.
//

#ifndef CGTOOL_GROOUTPUT_H
#define CGTOOL_GROOUTPUT_H

#include "trj_output.h"

class GROOutput : public TrjOutput{
protected:
    /** \brief Output file handle. */
    FILE *file_ = nullptr;

    /** \brief Open and prepare output file. */
    int openFile(const std::string &filename);
    /** \brief Close output file. */
    int closeFile();
public:
    /** \brief Constructor.  Calls openFile() .*/
    GROOutput(const int natoms, const std::string &filename);
    /** \brief Destructor.  Calls closeFile(). */
    ~GROOutput();

    /** \brief Write a Frame to output file. */
    int writeFrame(const Frame &frame);

    friend class Frame;
};

#endif //CGTOOL_GROOUTPUT_H
