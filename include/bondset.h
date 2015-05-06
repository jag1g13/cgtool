#ifndef BONDSET_H_
#define BONDSET_H_

#include <vector>
#include <string>
#include <map>

#include "bond_struct.h"
#include "frame.h"

using std::vector;
using std::string;

/**
* \brief Class that holds all bond lengths, angles and dihedrals to be calculated
*/
class BondSet{
protected:
    /** How many frames did we successfully measure */
    int numMeasures_ = 0;
    /** How many residues are being mapped using the same mapping? */
    int numResidues_ = 1;
    /** Temperature of simulation - for Boltzmann Inversion */
    double temp_ = 310.;
    /** Map of bead names to number - to put numbers into BondStructs */
    std::map<std::string, int> beadNums_;
    /** Number of threads - used for timing */
    int numThreads_ = 1;
    /** Residue name - for CSV output */
    std::string resname_;
    /** The residue that's present - to be made plural soon */
    Res residue_;

public:
    /** Vector of bond length pairs; Contains all bond lengths that must be calculated */
    vector<BondStruct> bonds_;
    /** Vector of bond angle triples */
    vector<BondStruct> angles_;
    /** Vector of bond dihedral quads */
    vector<BondStruct> dihedrals_;

    /** \brief Blank constructor */
    BondSet(){};

    /** \brief Constructor to read from file */
    BondSet(const std::string &cfgname, const int num_threads=1);

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

    /** \brief Perform Boltzmann Inversion on all bond_structs. */
    void BoltzmannInversion();

    /** \brief Calculate bond averages without full Boltzmann Inversion */
    void calcAvgs();

    /** \brief Write all bond parameters to CSVs.
    * SLOW.  This takes about the same amount of time as the complete
    * XTC input -> Boltzmann Inversion process so is turned off by default. */
    void writeCSV();
};

#endif