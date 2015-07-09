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

BondSet::BondSet(const string &cfgname, const vector<Residue> *residues,
                 const PotentialType potentials[3]){
    residues_ = residues;
    fromFile(cfgname);
}

void BondSet::fromFile(const string &filename){
    vector<string> tokens;
    Parser parser(filename);

    if(parser.getLineFromSection("temp", tokens, 1)) temp_ = stof(tokens[0]);

    if(parser.findSection("mapping")){
        int i = 0;
        while(parser.getLineFromSection("mapping", tokens, 1)){
            beadNums_[tokens[0]] = i;
            i++;
        }
    }else{
        beadNums_ = (*residues_)[0].name_to_num;
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
    for(int i=0; i < (*residues_)[0].num_residues; i++){
        bool res_okay = true;
        const int offset = i * (*residues_)[0].num_atoms;
        // Does the structure cross a pbc - will break bond lengths
        for(BondStruct &bond : bonds_){
            double dist = bond.bondLength(frame, offset);
            // If any distances > 1nm, molecule is on PBC
            // Also check for errors
            if(dist > 1. || std::isinf(dist) || std::isnan(dist)){
                res_okay = false;
                break;
            }
        }
        // Molecule is broken, skip
        if(!res_okay) continue;

        for(BondStruct &bond : bonds_){
            const double val = bond.bondLength(frame, offset);
            if(!std::isinf(val) && !std::isnan(val)) bond.values_.push_back(val);
        }
        for(BondStruct &bond : angles_){
            const double val = bond.bondAngle(frame, offset);
            if(!std::isinf(val) && !std::isnan(val)) bond.values_.push_back(val);
        }
        for(BondStruct &bond : dihedrals_){
            const double val = bond.bondAngle(frame, offset);
            if(!std::isinf(val) && !std::isnan(val)) bond.values_.push_back(val);
        }
        numMeasures_++;
    }
}

// Angles can't just be averaged like this - they wrap around
// Fine as approximation though, we won't deal much with angles close to 0
void BondSet::BoltzmannInversion(){
    if(numMeasures_ > 0){
        printf("Measured %'d molecules\n", numMeasures_);
    }else{
        printf("No bonds measured\n");
        return;
    }
    BoltzmannInverter bi(temp_);
    for(BondStruct &bond : bonds_) bi.calculate(bond);
    for(BondStruct &bond : angles_) bi.calculate(bond);
    for(BondStruct &bond : dihedrals_) bi.calculate(bond);
}

void BondSet::calcAvgs(){
    if(numMeasures_ > 0){
        printf("Measured %'d molecules\n", numMeasures_);
    }else{
        printf("No bonds measured\n");
        return;
    }
    BoltzmannInverter bi(temp_);
    for(BondStruct &bond : bonds_) bond.avg_ = bi.statisticalMoments(bond.values_);
    for(BondStruct &bond : angles_) bond.avg_ = bi.statisticalMoments(bond.values_);
    for(BondStruct &bond : dihedrals_) bond.avg_ = bi.statisticalMoments(bond.values_);
}

void BondSet::writeCSV(const int num_molecules) const{
    const string bond_file = (*residues_)[0].resname + "_bonds.dat";
    const string angle_file = (*residues_)[0].resname + "_angles.dat";
    const string dihedral_file = (*residues_)[0].resname + "_dihedrals.dat";

    backup_old_file(bond_file);
    backup_old_file(angle_file);
    backup_old_file(dihedral_file);

    FILE *f_bond = fopen(bond_file.c_str(), "w");
    FILE *f_angle = fopen(angle_file.c_str(), "w");
    FILE *f_dihedral = fopen(dihedral_file.c_str(), "w");

    // Scale increment so that ~num_molecules molecules are printed to CSV
    // Should be enough to be a good sample - but is much quicker than printing all
    int scale = 1;
    if(num_molecules > 0 && numMeasures_ > num_molecules)
        scale = static_cast<int>(numMeasures_ / static_cast<double>(num_molecules));

    for(int i=0; i < numMeasures_; i+=scale){
        for(const BondStruct &bond : bonds_) fprintf(f_bond, "%8.3f", bond.values_[i]);
        fprintf(f_bond, "\n");

        for(const BondStruct &bond : angles_) fprintf(f_angle, "%8.3f", bond.values_[i]);
        fprintf(f_angle, "\n");

        for(const BondStruct &bond : dihedrals_) fprintf(f_dihedral, "%8.3f", bond.values_[i]);
        fprintf(f_dihedral, "\n");
    }
    printf("Written %'d molecules to CSV\n", numMeasures_/scale);

    fclose(f_bond);
    fclose(f_angle);
    fclose(f_dihedral);
}