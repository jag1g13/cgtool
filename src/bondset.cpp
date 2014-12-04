#include "bondset.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include <boost/algorithm/string.hpp>

#include "parser.h"

#define DEBUG false

using std::vector;
using std::string;
using std::cout;
using std::endl;

BondStruct::BondStruct(int size){
    atom_names.resize(size);
    atom_nums.resize(size);
}

BondSet::BondSet(){
}

void BondSet::fromFile(string filename){
    vector<string> substrs;
    string section;
    Parser parser(filename);
    while(parser.getLine(&section, &substrs)){
        BondStruct bond_tmp = BondStruct(4);
        bond_tmp.atom_names = substrs;
        if(section == "length") {
            bonds_.push_back(bond_tmp);
        }else if(section == "angle") {
            angles_.push_back(bond_tmp);
        }else if(section == "dihedral"){
            dihedrals_.push_back(bond_tmp);
        }else{
            cout << "Error in bond sections" << endl;
        }
    }
    if(DEBUG){
        for(auto &i : dihedrals_){
            for(auto &j : i.atom_names){
                std::cout << " " << j;
            }
            std::cout << std::endl;
        }
    }
}

vector<float> BondSet::calcBondLens(Frame *frame){
    vector<float> bond_lens;
    for(vector<BondStruct>::iterator bond = bonds_.begin(); bond != bonds_.end(); ++bond) {
        bond_lens.push_back(frame->bondLength(&*bond));
    }
    return bond_lens;
}

vector<float> BondSet::calcBondAngles(Frame *frame){
    vector<float> bond_angles;
    for(vector<BondStruct>::iterator bond = angles_.begin(); bond != angles_.end(); ++bond){
        bond_angles.push_back(frame->bondAngle(&*bond));
    }
    return bond_angles;
}

vector<float> BondSet::calcBondDihedrals(Frame *frame){
    vector<float> bonds;
    for(vector<BondStruct>::iterator bond = dihedrals_.begin(); bond != dihedrals_.end(); ++bond){
        bonds.push_back(frame->bondAngle(&*bond));
    }
    return bonds;
}
