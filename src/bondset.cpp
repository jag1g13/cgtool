#include "bondset.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include <boost/algorithm/string.hpp>

#include "parser.h"
#include "frame.h"

using std::vector;
using std::string;
using std::cout;
using std::endl;

void BondSet::fromFile(string filename){
    vector<string> substrs;
    Parser parser(filename);
    while(parser.getLineFromSection("length", &substrs)) {
        BondStruct bond_tmp = BondStruct(4);
        bond_tmp.atom_names = substrs;
        bonds_.push_back(bond_tmp);
    }
    while(parser.getLineFromSection("angle", &substrs)) {
        BondStruct bond_tmp = BondStruct(4);
        bond_tmp.atom_names = substrs;
        angles_.push_back(bond_tmp);
    }
    while(parser.getLineFromSection("dihedral", &substrs)) {
        BondStruct bond_tmp = BondStruct(4);
        bond_tmp.atom_names = substrs;
        dihedrals_.push_back(bond_tmp);
    }
}

void BondSet::calcBondsInternal(Frame &frame){
    for(BondStruct &bond : bonds_){
        // does the structure cross a pbc - will break bond lengths
        float dist = frame.bondLength(bond);
        if(dist > 0.8f * frame.box_[0][0]){
            frame.invalid_ = true;
            return;
        }
        bond.values.push_back(dist);
    }
    for(BondStruct &bond : angles_){
        bond.values.push_back(frame.bondAngle(bond));
    }
    for(BondStruct &bond : dihedrals_){
        bond.values.push_back(frame.bondAngle(bond));
    }
}

void BondSet::calcAvgs(){
    for(BondStruct &bond : bonds_){
        for(float &b : bond.values){
            bond.avg += b;
        }
        bond.avg /= bond.values.size();
    }
    for(BondStruct &bond : angles_){
        for(float &b : bond.values){
            bond.avg += b;
        }
        bond.avg /= bond.values.size();
    }
    for(BondStruct &bond : dihedrals_){
        for(float &b : bond.values){
            bond.avg += b;
        }
        bond.avg /= bond.values.size();
    }
}
