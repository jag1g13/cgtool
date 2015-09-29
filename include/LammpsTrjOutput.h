//
// Created by james on 27/08/15.
//

#ifndef CGTOOL_LAMMPSTRJOUTPUT_H
#define CGTOOL_LAMMPSTRJOUTPUT_H

#include "trj_output.h"

class LammpsTrjOutput : public TrjOutput{
protected:
    /** \brief Output file handle. */
    FILE *file_ = nullptr;

    /** \brief Open and prepare output file. */
    int openFile(const std::string &filename);
    /** \brief Close output file. */
    int closeFile();
public:
    /** \brief Constructor.  Calls openFile() .*/
    LammpsTrjOutput(const int natoms, const std::string &filename);
    /** \brief Destructor.  Calls closeFile(). */
    ~LammpsTrjOutput();

    /** \brief Write a Frame to output file. */
    int writeFrame(const Frame &frame);

    friend class Frame;
};


#endif //CGTOOL_LAMMPSTRJOUTPUT_H
