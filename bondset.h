
using std::vector;
using std::string;

#ifndef BONDSTRUCT
#define BOND_STRUCT
/**
* \brief Struct to hold atoms in bonds, angles and dihedrals
*/
struct bond_struct{
    vector<string> atom_names;
    vector<int> atom_nums;
};
#endif

/**
* \brief Class that holds all bond lengths, angles and dihedrals to be calculated
*/
class BondSet{
public:
    vector<bond_struct> bonds, angles, dihedrals;

    BondSet();

    bool from_file(string);
};