#ifndef BONDSET_H
#define BONDSET_H

#include <vector>
#include <string>

using std::vector;
using std::string;

/**
* \brief Struct to hold atoms in bonds, angles and dihedrals
*/
struct BondStruct{
    /** Vector of atom names for this bond property; For a bond length will contain two names; three for angle; four for dihedral */
    vector<string> atom_names;
    /** Vector of atom numbers for this bond property; For a bond length will contain two names; three for angle; four for dihedral */
    vector<int> atom_nums;
    /** Constructor to set size (bond/angle/dihedral) */
    BondStruct(int size);
};

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
    bool fromFile(string);
};

#endif