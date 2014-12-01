#ifndef CGMAP_H
#define CGMAP_H

#include <string>
#include <vector>

#include "frame.h"

using std::string;
using std::vector;

/**
* \brief Struct to hold the mapping for a single CG bead
*/
struct BeadMap{
    /** The name of this CG bead */
    string cg_bead;
    /** The atoms which should be mapped into this bead */
    vector<string> atoms;
};

/**
* \brief Contains data and functions related to the CG mapping
*
* Has functions to read in a CG mapping from file and apply it to an atomistic Frame
* Mostly just a wrapper around a BeadMap vector
*/
class CGMap{
public:
    /** Vector of BeadMap; holds the mappings for every bead */
    vector<BeadMap> mapping;
    /** Number of beads defined */
    int num_beads;

    /**
    * \brief Constructor to create a blank instance
    */
    CGMap();

    /**
    * \brief Constructor to create an instance from the mapping file provided
    */
    CGMap(string filename);

    /**
    * \brief Read in CG mapping from file
    */
    bool from_file(string filename);

    /**
    * \brief Setup a CG Frame object that has already been declared
    */
    void init_frame(Frame *cg_frame);

    /**
    * \brief Apply CG mapping to an atomistic Frame
    *
    * Requires a pre-constructed Frame for output, but can allocate the number of atoms here
    */
    bool apply(Frame *aa_frame, Frame *cg_frame);
};

#endif