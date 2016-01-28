#ifndef BONDSET_H_
#define BONDSET_H_

#include <vector>
#include <string>
#include <map>

#include "bond_struct.h"
#include "frame.h"
#include "file_io.h"

using std::vector;
using std::string;

/**
* \brief Class that holds all bond lengths, angles and dihedrals to be calculated
*/
class BondSet{
protected:
    /** How many frames did we successfully measure */
    int numMeasures_ = 0;
    /** Temperature of simulation - for Boltzmann Inversion */
    double temp_ = 310.;
    /** Map of bead names to number */
    std::map<std::string, int> beadNums_;
    /** The residues */
    const std::vector<Residue> &residues_;

public:
    /** Vector of bond length pairs; Contains all bond lengths that must be calculated */
    vector<BondStruct> bonds_;
    /** Vector of bond angle triples */
    vector<BondStruct> angles_;
    /** Vector of bond dihedral quads */
    vector<BondStruct> dihedrals_;

    /** \brief Constructor to read from file */
    BondSet(const std::string &cfgname, const std::vector<Residue> &residues,
            const PotentialType potentials[3], const double temp);

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
    void writeCSV(const int num_molecules) const;
};

#endif