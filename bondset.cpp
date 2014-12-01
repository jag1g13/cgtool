#include "bondset.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <boost/algorithm/string.hpp>

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

bool BondSet::fromFile(string filename){
    bool ok = true;
    std::ifstream map_file(filename);
    string line;
    vector<string> substrs;
    if(!map_file.is_open()) return 0;   // couldn't open file
    while(getline(map_file, line)){
        if(line[0] == ';' || line[0] == '#') continue;  // skip comments
        if(line[0] == '[') continue;                    // TODO change this
        if(line == "") continue;                        // line is empty, ignore it
        boost::split(substrs, line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
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
    return ok;
}
