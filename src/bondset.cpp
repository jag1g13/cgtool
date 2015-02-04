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
using std::fprintf;

void BondSet::fromFile(const string &filename){
    //TODO don't read in anything if there isn't a line - why does it do this?
    vector<string> substrs;
    Parser parser(filename);
    while(parser.getLineFromSection("length", substrs)) {
        BondStruct bond_tmp = BondStruct(4);
        bond_tmp.atomNames_ = substrs;
        bonds_.push_back(bond_tmp);
    }
    while(parser.getLineFromSection("angle", substrs)) {
        BondStruct bond_tmp = BondStruct(4);
        bond_tmp.atomNames_ = substrs;
        angles_.push_back(bond_tmp);
    }
    while(parser.getLineFromSection("dihedral", substrs)) {
        BondStruct bond_tmp = BondStruct(4);
        bond_tmp.atomNames_ = substrs;
        dihedrals_.push_back(bond_tmp);
    }
}

void BondSet::calcBondsInternal(Frame &frame){
    for(BondStruct &bond : bonds_){
        // does the structure cross a pbc - will break bond lengths
        float dist = frame.bondLength(bond);
//        if(dist > 0.8f * frame.box_[0][0]){
        if(dist > 1.f || dist < 0.01f){
            frame.invalid_ = true;
            return;
        }
    }
    numFrames_++;
    for(BondStruct &bond : bonds_){
        bond.values_.push_back(frame.bondLength(bond));
    }
    for(BondStruct &bond : angles_){
        bond.values_.push_back(frame.bondAngle(bond));
    }
    for(BondStruct &bond : dihedrals_){
        bond.values_.push_back(frame.bondAngle(bond));
    }
}

void BondSet::calcAvgs(){
    for(BondStruct &bond : bonds_){
        bond.calcAvg();
        bond.binHistogram(100);
    }
    for(BondStruct &bond : angles_){
        bond.calcAvg();
    }
    for(BondStruct &bond : dihedrals_){
        bond.calcAvg();
    }
}

void BondSet::writeCSV(){
    FILE *f_bond = fopen("bonds.csv", "w");
    FILE *f_angle = fopen("angles.csv", "w");
    FILE *f_dihedral = fopen("dihedrals.csv", "w");
    for(int i=0; i < numFrames_; i++){
        for(BondStruct &bond : bonds_){
            fprintf(f_bond, "%8.4f", bond.values_[i]);
        }
        fprintf(f_bond, "\n");
        for(BondStruct &bond : angles_){
            fprintf(f_angle, "%8.4f", bond.values_[i]);
        }
        fprintf(f_angle, "\n");
        for(BondStruct &bond : dihedrals_){
            fprintf(f_dihedral, "%8.4f", bond.values_[i]);
        }
        fprintf(f_dihedral, "\n");
    }
    fclose(f_bond);
    fclose(f_angle);
    fclose(f_dihedral);
}