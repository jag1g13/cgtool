#ifndef BONDSTRUCT_H_
#define BONDSTRUCT_H_

#include <vector>
#include <string>

enum class BondType{LENGTH, ANGLE, DIHEDRAL};

/**
* \brief Class to hold atoms in bonds, angles and dihedrals.
*/
class BondStruct{
public:
    /** The values of the bond parameter (length, angle, dih) for each Frame */
    std::vector<double> values_;
    /** Vector of atom numbers for this bond property; For a bond length will contain two names; three for angle; four for dihedral */
    std::vector<int> atomNums_;
    /** Average of bond parameter.  Double avoids overflow on summing over large numbers of Frames */
    double avg_ = 0.;
    /** Force constant */
    double forceConstant_ = 0.;
    /** \brief R^2 of fitting gaussian to bond distribution
    * A low value indicates that the bond is probably bimodal */
    double rsqr_ = 0.;
    /** What type of bond measure is it?  Length, angle or dihedral */
    BondType type_;

    /** Constructor to set size (bond/angle/dihedral) */
    BondStruct(const int size);
    /** Blank constructor */
    BondStruct(){};
    /** Copy constructor - required to push into vectors */
    BondStruct(const BondStruct &other);
};

#endif