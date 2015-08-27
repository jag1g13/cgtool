//
// Created by james on 27/08/15.
//

#ifndef CGTOOL_LAMMPSDATAOUTPUT_H
#define CGTOOL_LAMMPSDATAOUTPUT_H

#include "trj_output.h"

class LammpsDataOutput : public TrjOutput{
protected:
    /** \brief Output file handle. */
    FILE *file_ = nullptr;

    /** \brief Open and prepare output file. */
    int openFile(const std::string &filename);
    /** \brief Close output file. */
    int closeFile();
public:
    /** \brief Constructor.  Calls openFile() .*/
    LammpsDataOutput(const int natoms, const std::string &filename);
    /** \brief Destructor.  Calls closeFile(). */
    ~LammpsDataOutput();

    /** \brief Write a Frame to output file. */
    int writeFrame(const Frame &frame);

    friend class Frame;
};


#endif //CGTOOL_LAMMPSDATAOUTPUT_H
