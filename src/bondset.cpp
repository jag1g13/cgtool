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
    vector<string> substrs;
    Parser parser(filename);

    int i = 0;
    while(parser.getLineFromSection("mapping", substrs)){
//        beadNums_.insert({substrs[0], i});
        beadNums_.emplace(substrs[0], i);
        cout << substrs[0] << " " << beadNums_[substrs[0]] << endl;
    }

    while(parser.getLineFromSection("length", substrs)){
        bonds_.emplace_back(BondStruct(2));
        for(int j = 0; j < 2; j++){
            bonds_.back().atomNames_[j] = substrs[j];
            //TODO these numbers don't work - all are 0
            bonds_.back().atomNums_[j] = beadNums_[substrs[j]];
            cout << substrs[j] << " " << beadNums_[substrs[j]] << endl;
        }
//        bonds_.back().atomNames_ = substrs;
//        bonds_.back().atomNums_.emplace_back(beadNums_[name]);
    }
    while(parser.getLineFromSection("angle", substrs)){
        angles_.emplace_back(BondStruct(3));
        for(int j = 0; j < 3; j++){
            angles_.back().atomNames_[j] = substrs[j];
            angles_.back().atomNums_[j] = beadNums_[substrs[j]];
        }
//        angles_.back().atomNames_ = substrs;
//        for(const string &name : angles_.back().atomNames_){
//            angles_.back().atomNums_.emplace_back(beadNums_[name]);
//        }
    }
    while(parser.getLineFromSection("dihedral", substrs)){
        dihedrals_.emplace_back(BondStruct(4));
        for(int j = 0; j < 4; j++){
            dihedrals_.back().atomNames_[j] = substrs[j];
            dihedrals_.back().atomNums_[j] = beadNums_[substrs[j]];
        }
//        dihedrals_.back().atomNames_ = substrs;
//        for(const string &name : dihedrals_.back().atomNames_){
//            dihedrals_.back().atomNums_.emplace_back(beadNums_[name]);
//        }
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