#ifndef BONDSET_H_
#define BONDSET_H_

#include <vector>
#include <string>

#include "frame.h"
#include "bond_struct.h"

using std::vector;
using std::string;


/**
* \brief Class that holds all bond lengths, angles and dihedrals to be calculated
*/
class BondSet{
public:
    /** Vector of bond length pairs; Contains all bond lengths that must be calculated */
    vector<BondStruct> bonds_;
    /** Vector of bond angle triples */
    vector<BondStruct> angles_;
    /** Vector of bond dihedral quads */
    vector<BondStruct> dihedrals_;

    BondSet();

    /**
    * \brief Reads in from file all bond properties to be calculated
    *
    * Gets Vectors of all bond lengths, angles and dihedrals that must be calculated.
    */
    void fromFile(string);

    /**
    * \brief Calculate all bond lengths that were requested in the input file
    *
    * *Something, something, locality*
    */
    vector<float> calcBondLens(Frame *frame);
};

#endif