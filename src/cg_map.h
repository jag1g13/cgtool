#ifndef CGMAP_H_
#define CGMAP_H_

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
    /** The number of atoms this bead contains */
    int num_atoms;
    /** The atoms which should be mapped into this bead */
    vector<string> atoms;
    /** The atoms which should be mapped into this bead, by order in Frame */
    vector<int> atom_nums;
    /** Total mass of bead */
    float mass;
    /** Total charge on bead */
    float charge;
};

enum class MapType{CM, GC, ATOM};


/**
* \brief Contains data and functions related to the CG mapping
*
* Has functions to read in a CG mapping from file and apply it to an atomistic Frame
* Mostly just a wrapper around a BeadMap vector
*/
class CGMap{
public:
    /** Vector of BeadMap; holds the mappings for every bead */
    vector<BeadMap> mapping_;
    /** Dictionary mapping an atom to the bead it should be mapped into */
    std::map<string, BeadMap*> atomname_to_bead_;
    /** Number of beads defined */
    int num_beads;
    /** \brief What type of mapping are we going to apply?  CM, GC, or atom centred
    * Default is mapping on first atom in bead. */
    MapType mapType_ = MapType::ATOM;

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
    *
    * Will throw std::runtime_error if file cannot be opened
    */
    void fromFile(string filename);

    /**
    * \brief Setup a CG Frame object that has already been declared
    *
    * Allocates space for each bead and copies over constant data from the atomistic Frame
    */
    Frame initFrame(const Frame &aa_frame);

    /**
    * \brief Apply CG mapping to an atomistic Frame
    *
    * Requires that initFrame has already been called to setup the CG Frame.
    * Will throw std::runtime_error if Frame hasn't been setup.
    */
    bool apply(const Frame &aa_frame, Frame &cg_frame);
};

#endif