#ifndef CGMAP_H
#define CGMAP_H

#include <string>

#include "frame.h"

using std::string;

/**
* \brief Contains data and functions related to the CG mapping
*
* Has functions to read in a CG mapping from file and apply it to an atomistic Frame
*/
class CGMap{
public:

    CGMap();

    /**
    * \brief Read in CG mapping from file
    */
    bool from_file(string filename);

    /**
    * \brief Apply CG mapping to an atomistic Frame
    *
    * Requires a pre-constructed Frame for output, but can allocate the number of atoms here
    */
    bool apply(Frame* aa_frame, Frame* cg_frame);
};

#endif