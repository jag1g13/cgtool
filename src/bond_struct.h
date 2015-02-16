#ifndef BONDSTRUCT_H_
#define BONDSTRUCT_H_

#include <vector>
#include <string>

#include "array.h"

/**
* \brief Class to hold atoms in bonds, angles and dihedrals.
*/
class BondStruct{
protected:
    /** Keep track of serial number */
    static int totalNum_;
    /** Serial number, count includes all bond lengths, angles and dihedrals */
    int num_;

public:
    /** The values of the bond parameter (length, angle, dih) for each Frame */
    std::vector<double> values_;
    /** Vector of atom numbers for this bond property; For a bond length will contain two names; three for angle; four for dihedral */
    std::vector<int> atomNums_;
    /** Average of bond parameter.  Double avoids overflow on summing over large numbers of Frames */
    double avg_ = 0.;

    /** Constructor to set size (bond/angle/dihedral) */
    BondStruct(const int size);
    /** Blank constructor */
    BondStruct(){num_ = totalNum_; totalNum_++;};
    /** Copy constructor - required to push into vectors */
    BondStruct(const BondStruct &other);

    /** Calculate average of the bonds */
    void calcAvg();
};

#endif
