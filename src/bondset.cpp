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

BondStruct::BondStruct(const BondStruct &other){
    unsigned long size = other.atomNames_.size();
    atomNames_.resize(size);
    atomNums_.resize(size);
    for(int i = 0; i < size; i++){
        atomNames_[i] = other.atomNames_[i];
        atomNums_[i] = other.atomNums_[i];
    }
}

void BondSet::fromFile(const string &filename){
    vector<string> substrs;
    Parser parser(filename);
    while(parser.getLineFromSection("length", substrs)) {
        bonds_.emplace_back(BondStruct(2));
        bonds_.back().atomNames_ = substrs;
    }
    while(parser.getLineFromSection("angle", substrs)) {
        angles_.emplace_back(BondStruct(3));
        angles_.back().atomNames_ = substrs;
    }
    while(parser.getLineFromSection("dihedral", substrs)) {
        dihedrals_.emplace_back(BondStruct(4));
        dihedrals_.back().atomNames_ = substrs;
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
//        bond.binHistogram(100);
    }
    for(BondStruct &bond : angles_){
        bond.calcAvg();
//        bond.binHistogram(100);
    }
    for(BondStruct &bond : dihedrals_){
        bond.calcAvg();
//        bond.binHistogram(100);
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