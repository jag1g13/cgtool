#ifndef BONDSTRUCT_H_
#define BONDSTRUCT_H_

#include <vector>
#include <string>

#include "array.h"

/**
* \brief Class to hold atoms in bonds, angles and dihedrals.
* Can calculate stats and perform Boltzmann Inversion.
*/
class BondStruct{
public:
    /** Vector of atom names for this bond property; For a bond length will contain two names; three for angle; four for dihedral */
    std::vector<std::string> atomNames_;
    /** Vector of atom numbers for this bond property; For a bond length will contain two names; three for angle; four for dihedral */
    std::vector<int> atomNums_;
    /** The values of the bond parameter (length, angle, dih) for each Frame */
    std::vector<float> values_;
    /** Average bond parameter.  Double stops overflow on summing */
    double avg_ = 0.;

    /** Constructor to set size (bond/angle/dihedral) */
    BondStruct(const int size);
    /** Blank constructor */
    BondStruct(){};
    /** Copy constructor - required to push into vectors */
    BondStruct(const BondStruct &other);

    /** Calculate average of the bonds */
    void calcAvg();
};

#endif
