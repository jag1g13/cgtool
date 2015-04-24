#include "bondset.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <ctime>
#include <cmath>

#include <boost/algorithm/string.hpp>

#include "parser.h"
#include "boltzmann_inverter.h"
#include "small_functions.h"

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::fprintf;

BondSet::BondSet(const std::string &cfgname, const int num_threads){
    fromFile(cfgname);
    numThreads_ = num_threads;
}

void BondSet::fromFile(const string &filename){
    vector<string> tokens;
    Parser parser(filename);

    if(parser.getLineFromSection("residues", tokens, 2)) numResidues_ = stoi(tokens[0]);
    if(parser.getLineFromSection("temp", tokens, 1)) temp_ = stof(tokens[0]);

    int i = 0;
    while(parser.getLineFromSection("mapping", tokens, 3)){
        beadNums_.emplace(tokens[0], i++);
    }

    //TODO Can emplace_back() be replaced?
    while(parser.getLineFromSection("length", tokens, 2)){
        bonds_.emplace_back(BondStruct(2));
        bonds_.back().atomNums_[0] = beadNums_[tokens[0]];
        bonds_.back().atomNums_[1] = beadNums_[tokens[1]];
    }

    while(parser.getLineFromSection("angle", tokens, 3)){
        angles_.emplace_back(BondStruct(3));
        angles_.back().atomNums_[0] = beadNums_[tokens[0]];
        angles_.back().atomNums_[1] = beadNums_[tokens[1]];
        angles_.back().atomNums_[2] = beadNums_[tokens[2]];
    }

    while(parser.getLineFromSection("dihedral", tokens, 4)){
        dihedrals_.emplace_back(BondStruct(4));
        dihedrals_.back().atomNums_[0] = beadNums_[tokens[0]];
        dihedrals_.back().atomNums_[1] = beadNums_[tokens[1]];
        dihedrals_.back().atomNums_[2] = beadNums_[tokens[2]];
        dihedrals_.back().atomNums_[3] = beadNums_[tokens[3]];
    }
}

void BondSet::calcBondsInternal(Frame &frame){
    for(int i=0; i < numResidues_; i++){
        bool res_okay = true;
        const int offset = i * frame.numAtomsPerResidue_;
        // Does the structure cross a pbc - will break bond lengths
        for(BondStruct &bond : bonds_){
            double dist = frame.bondLength(bond, offset);
            // If any distances > 1nm, molecule is on PBC
            if(dist > 2.f || dist < 0.01f){
                res_okay = false;
                break;
            }
        }
        // Molecule is broken, skip
        if(!res_okay) continue;

        for(BondStruct &bond : bonds_){
            const double val = frame.bondLength(bond, offset);
            if(!std::isinf(val)) bond.values_.push_back(val);
        }
        for(BondStruct &bond : angles_){
            const double val = frame.bondAngle(bond, offset);
            if(!std::isinf(val)) bond.values_.push_back(val);
        }
        for(BondStruct &bond : dihedrals_){
            const double val = frame.bondAngle(bond, offset);
            if(!std::isinf(val)) bond.values_.push_back(val);
        }
        numMeasures_++;
    }
}

void BondSet::BoltzmannInversion(){
    if(numMeasures_ > 0){
        cout << "Measured " << numMeasures_ << " molecules" << endl;
    }else{
        cout << "No bonds measured" << endl;
        return;
    }
    BoltzmannInverter bi(temp_);
    for(BondStruct &bond : bonds_) bi.calculate(bond);
    for(BondStruct &bond : angles_) bi.calculate(bond);
    for(BondStruct &bond : dihedrals_) bi.calculate(bond);
}

void BondSet::calcAvgs(){
    if(numMeasures_ > 0){
        cout << "Measured " << numMeasures_ << " molecules" << endl;
    }else{
        cout << "No bonds measured" << endl;
        return;
    }
    BoltzmannInverter bi(temp_);
    for(BondStruct &bond : bonds_) bond.avg_ = bi.statisticalMoments(bond.values_);
    for(BondStruct &bond : angles_) bond.avg_ = bi.statisticalMoments(bond.values_);
    for(BondStruct &bond : dihedrals_) bond.avg_ = bi.statisticalMoments(bond.values_);
}

void BondSet::writeCSV(){
    FILE *f_bond = fopen("bonds.csv", "w");
    FILE *f_angle = fopen("angles.csv", "w");
    FILE *f_dihedral = fopen("dihedrals.csv", "w");
    clock_t start = std::clock();

    // Scale increment so that ~10k molecules are printed to CSV
    // Should be enough to be a good sample - but is much quicker
    int scale = 1;
    if(numMeasures_ > 1e4) scale = numMeasures_ / 1e4;

    for(int i=0; i < numMeasures_; i+=scale){
        #ifdef UPDATE_PROGRESS
        if(i % 50000 == 0){
            float time = time_since(start, numThreads_);
            float fps = i / time;
            float t_remain = (numMeasures_ - i) / fps;
            printf("Written %d molecules to CSV %6.1fs remaining\r", i, t_remain);
            std::flush(cout);
        }
        #endif

        for(BondStruct &bond : bonds_) fprintf(f_bond, "%8.3f", bond.values_[i]);
        fprintf(f_bond, "\n");

        for(BondStruct &bond : angles_) fprintf(f_angle, "%8.3f", bond.values_[i]);
        fprintf(f_angle, "\n");

        for(BondStruct &bond : dihedrals_) fprintf(f_dihedral, "%8.3f", bond.values_[i]);
        fprintf(f_dihedral, "\n");
    }
    cout << string(80, ' ') << "\r";
    cout << "Written  " << numMeasures_/scale << " molecules to CSV" << endl;

    fclose(f_bond);
    fclose(f_angle);
    fclose(f_dihedral);
}