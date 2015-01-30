#ifndef BOND_STRUC_H_
#define BOND_STRUC_H_

#include <vector>
#include <string>

/**
* \brief Struct to hold atoms in bonds, angles and dihedrals
*/
struct BondStruct{
    /** Vector of atom names for this bond property; For a bond length will contain two names; three for angle; four for dihedral */
    std::vector<std::string> atom_names;
    /** Vector of atom numbers for this bond property; For a bond length will contain two names; three for angle; four for dihedral */
    std::vector<int> atom_nums;
    /** The value of the bond parameter */
    float value;
    /** Constructor to set size (bond/angle/dihedral) */
    BondStruct(int size);
};

#endif
