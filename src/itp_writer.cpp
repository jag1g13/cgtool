#include "itp_writer.h"

#include <iostream>

using std::string;
using std::fprintf;
using std::cout;
using std::endl;
using std::vector;

ITPWriter::ITPWriter(string name){
    name_ = name;
    itp_ = std::fopen(name.c_str(), "w");
    if(itp_ == NULL){
        cout << "Could not open itp file for writing" << endl;
        exit(-1);
    }
    fprintf(itp_, "%s", header_.c_str());
}

void ITPWriter::newSection(std::string section_name){
    section_ = section_name;
    fprintf(itp_, "\n[ %s ]\n", section_.c_str());
}

void ITPWriter::printAtoms(const CGMap &map){
    newSection("atoms");
    for(BeadMap bead : map.mapping_){
        fprintf(itp_, "%6i %10s %6i %6s %6s %6i %10.4f %10.4f;\n",
                bead.num, bead.type.c_str(), bead.num, bead.name.c_str(),
                bead.name.c_str(), bead.num, bead.charge, bead.mass);
    }
}

//TODO use atom numbers instead of names to make a valid GROMACS topology/itp
void ITPWriter::printBonds(const BondSet &bond_set){
    newSection("bonds");
    for(BondStruct bond : bond_set.bonds_){
        // I don't know what the last integer is
//        fprintf(itp_, "%5i %5i %5i %12.5f %12.5e;\n",
//                bond.atom_nums[0], bond.atom_nums[1], 2,
//                bond.avg, 0.f);
        fprintf(itp_, "%5s %5s %5i %12.5f %12.5e;\n",
                bond.atom_names[0].c_str(), bond.atom_names[1].c_str(), 2,
                bond.avg, 0.f);
    }
    newSection("bonds");
    for(BondStruct bond : bond_set.angles_){
//        fprintf(itp_, "%5i %5i %5i %5i %12.5f %12.5e;\n",
//                bond.atom_nums[0], bond.atom_nums[1],
//                bond.atom_nums[2], 2,
//                bond.avg, 0.f);
        fprintf(itp_, "%5s %5s %5s %5i %12.5f %12.5e;\n",
                bond.atom_names[0].c_str(), bond.atom_names[1].c_str(),
                bond.atom_names[2].c_str(), 2,
                bond.avg, 0.f);
    }
    newSection("bonds");
    for(BondStruct bond : bond_set.dihedrals_){
//        fprintf(itp_, "%5i %5i %5i %5i %5i %12.5f %12.5e;\n",
//                bond.atom_nums[0], bond.atom_nums[1],
//                bond.atom_nums[2], bond.atom_nums[3],
//                2, bond.avg, 0.f);
        fprintf(itp_, "%5s %5s %5s %5s %5i %12.5f %12.5e;\n",
                bond.atom_names[0].c_str(), bond.atom_names[1].c_str(),
                bond.atom_names[2].c_str(), bond.atom_names[3].c_str(),
                2, bond.avg, 0.f);
    }
}