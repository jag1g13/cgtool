#include "itp_writer.h"

#include <iostream>

using std::string;
using std::fprintf;
using std::cout;
using std::endl;
using std::vector;

ITPWriter::ITPWriter(const string &name){
    name_ = name;
    itp_ = std::fopen(name.c_str(), "w");
    if(itp_ == NULL){
        cout << "Could not open itp file for writing" << endl;
        exit(-1);
    }
    fprintf(itp_, "%s", header_.c_str());
}

void ITPWriter::newSection(const string &section_name){
    section_ = section_name;
    fprintf(itp_, "\n[ %s ]\n", section_.c_str());
}

void ITPWriter::printAtoms(const CGMap &map){
    newSection("atoms");
    for(BeadMap bead : map.mapping_){
        fprintf(itp_, "%6i %10s %6i %6s %6s %6i %10.4f %10.4f;\n",
                bead.num+1, bead.type.c_str(), 1, bead.name.c_str(),
                bead.name.c_str(), bead.num+1, bead.charge, bead.mass);
    }
}

void ITPWriter::printBonds(const BondSet &bond_set){
    newSection("bonds");
    for(const BondStruct &bond : bond_set.bonds_){
        // I don't know what the last integer is
        fprintf(itp_, "%5i %5i %5i %12.5f %12.5e;\n",
                bond.atomNums_[0]+1, bond.atomNums_[1]+1, 1,
                bond.avg_, bond.forceConstant_);
    }
    newSection("angles");
    for(const BondStruct &bond : bond_set.angles_){
        fprintf(itp_, "%5i %5i %5i %5i %12.5f %12.5e;\n",
                bond.atomNums_[0]+1, bond.atomNums_[1]+1,
                bond.atomNums_[2]+1, 2,
                bond.avg_, bond.forceConstant_);
    }
    newSection("dihedrals");
    for(const BondStruct &bond : bond_set.dihedrals_){
        fprintf(itp_, "%5i %5i %5i %5i %5i %12.5f %12.5e;\n",
                bond.atomNums_[0]+1, bond.atomNums_[1]+1,
                bond.atomNums_[2]+1, bond.atomNums_[3]+1,
                1, bond.avg_, bond.forceConstant_);
    }
}