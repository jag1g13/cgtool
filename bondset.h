#ifndef BONDSET_H
#define BONDSET_H
using std::vector;
using std::string;

/**
* \brief Struct to hold atoms in bonds, angles and dihedrals
*/
struct bond_struct{
    vector<string> atom_names;
    vector<int> atom_nums;
};

/**
* \brief Class that holds all bond lengths, angles and dihedrals to be calculated
*/
class BondSet{
public:
    vector<bond_struct> bonds, angles, dihedrals;

    BondSet();

    bool from_file(string);
};

#endif