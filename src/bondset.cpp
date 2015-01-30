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
    string section;
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

vector<float> BondSet::calcBondLens(Frame &frame){
    vector<float> bonds;
    vector<float> empty;
    for(auto &bond : bonds_){
        bonds.push_back(frame.bondLength(&bond));
        // does the structure cross a pbc - will break bond lengths
//        if(*bonds.end() > 0.8f * frame->box_[0][0]){
//            frame->invalid_ = true;
//            return empty;
//        }
    }
    return bonds;
}

vector<float> BondSet::calcBondAngles(Frame &frame){
    vector<float> bonds;
    vector<float> empty;
    if(frame.invalid_) return empty;
    for(auto &bond : angles_){
        bonds.push_back(frame.bondAngle(&bond));
    }
    return bonds;
}

vector<float> BondSet::calcBondDihedrals(Frame &frame){
    vector<float> bonds;
    vector<float> empty;
    if(frame.invalid_) return empty;
    for(auto &bond : dihedrals_){
        bonds.push_back(frame.bondAngle(&bond));
    }
    return bonds;
}

void BondSet::calcBondsInternal(Frame &frame){
    for(BondStruct &bond : bonds_){
        // does the structure cross a pbc - will break bond lengths
        float dist = frame.bondLength(&bond);
        if(dist > 0.8f * frame.box_[0][0]){
            frame.invalid_ = true;
            return;
        }
        bond.values.push_back(frame.bondLength(&bond));
    }
    for(BondStruct &bond : angles_){
        bond.values.push_back(frame.bondAngle(&bond));
    }
    for(BondStruct &bond : dihedrals_){
        bond.values.push_back(frame.bondAngle(&bond));
    }
}
