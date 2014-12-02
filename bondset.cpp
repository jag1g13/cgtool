#include "bondset.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include <boost/algorithm/string.hpp>

#include "parser.h"

#define DEBUG true

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
        BondStruct bond_tmp = BondStruct(2);
        bond_tmp.atom_names = substrs;
        switch(substrs.size()){
            case 2:
                bonds_.push_back(bond_tmp);
                break;
            case 3:
                angles_.push_back(bond_tmp);
                break;
            case 4:
                dihedrals_.push_back(bond_tmp);
                break;
            default: cout << "Case error in BondSet::fromFile" << endl;
        }
    }
    if(DEBUG){
        for(auto &i : bonds_){
            for(auto &j : i.atom_names){
                std::cout << " " << j;
            }
            std::cout << std::endl;
        }
    }
}

vector<float> BondSet::calcBondLens(Frame *frame){
    vector<float> bond_lens(bonds_.size());
    for(vector<BondStruct>::iterator bond = bonds_.begin(); bond != bonds_.end(); ++bond){
        //bond_lens.push_back(frame->bondLength(bond->atom_nums[0], bond->atom_nums[1]));
        bond_lens.push_back(frame->bondLength(&*bond));
    }
    //cout << bond_lens.back() << endl;
    return bond_lens;
}
