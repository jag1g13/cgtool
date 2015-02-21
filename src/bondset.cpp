#include "bondset.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include <boost/algorithm/string.hpp>

#include "parser.h"
#include "boltzmann_inverter.h"

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
        beadNums_.emplace(substrs[0], i);
        i++;
    }

    while(parser.getLineFromSection("length", substrs)){
        bonds_.emplace_back(BondStruct(2));
        for(int j = 0; j < 2; j++){
            bonds_.back().atomNums_[j] = beadNums_[substrs[j]];
        }
    }
    while(parser.getLineFromSection("angle", substrs)){
        angles_.emplace_back(BondStruct(3));
        for(int j = 0; j < 3; j++){
            angles_.back().atomNums_[j] = beadNums_[substrs[j]];
        }
    }
    while(parser.getLineFromSection("dihedral", substrs)){
        dihedrals_.emplace_back(BondStruct(4));
        for(int j = 0; j < 4; j++){
            dihedrals_.back().atomNums_[j] = beadNums_[substrs[j]];
        }
    }
}

void BondSet::calcBondsInternal(Frame &frame){
    for(BondStruct &bond : bonds_){
//        does the structure cross a pbc - will break bond lengths
        float dist = frame.bondLength(bond);
//        if(dist > 0.8f * frame.box_[0][0])
        if(dist > 1.f || dist < 0.01f){
            frame.invalid_ = true;
            cout << "Molecule must not be broken by PBC" << endl;
            cout << "Consider using trjconv -pbc nojump" << endl;
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
    }
    for(BondStruct &bond : angles_){
        bond.calcAvg();
    }
    for(BondStruct &bond : dihedrals_){
        bond.calcAvg();
    }
}

void BondSet::boltzmannInversion(){
    BoltzmannInverter bi;
    for(BondStruct &bond : bonds_){
        bond.calcAvg();
        bi.statisticalMoments(bond.values_);
        bi.binHistogram(bond, 32);
        bi.gaussianRSquared();
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