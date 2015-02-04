#ifndef BONDSET_H_
#define BONDSET_H_

#include <vector>
#include <string>

using std::vector;
using std::string;

// resolve circular dependency
class Frame;

/**
* \brief Struct to hold atoms in bonds, angles and dihedrals
*/
struct BondStruct{
    /** Vector of atom names for this bond property; For a bond length will contain two names; three for angle; four for dihedral */
    std::vector<std::string> atom_names;
    /** Vector of atom numbers for this bond property; For a bond length will contain two names; three for angle; four for dihedral */
    std::vector<int> atom_nums;
    /** The values of the bond parameter (length, angle, dih) for each Frame */
    std::vector<float> values;
    /** Average bond parameter */
    double avg;
    /** Constructor to set size (bond/angle/dihedral) */
    BondStruct(int size){atom_names.resize(size); atom_nums.resize(size);};
    BondStruct(){};
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
    /** How many frames did we successfully measure */
    int numFrames_ = 0;

    /** \brief Blank constructor */
    BondSet(){};

    /** \brief Constructor to read from file */
    BondSet(const std::string &cfgname){fromFile(cfgname);};

    /**
    * \brief Reads in from file all bond properties to be calculated
    *
    * Gets Vectors of all bond lengths, angles and dihedrals that must be calculated.
    */
    void fromFile(const string &filename);

    /**
    * \brief Calculate all bond lengths, angles and dihedrals.
    * There are stored inside the BondStructs to be passed to averaging functions later.
    */
    void calcBondsInternal(Frame &frame);

    /** \brief Calculate averages of all bond measurements */
    void calcAvgs();

    /** \brief Write all values to CSVs */
    void writeCSV();
};

#endif