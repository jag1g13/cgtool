#include "itp_writer.h"

#include <iostream>
#include <ctime>

#include "small_functions.h"

using std::string;
using std::fprintf;
using std::cout;
using std::endl;
using std::vector;

ITPWriter::ITPWriter(const string &resName, string itpname){
    if(itpname == "") itpname = resName + "CG.itp";
    name_ = itpname;
    resName_ = resName;
    backup_old_file(name_.c_str());
    itp_ = std::fopen(name_.c_str(), "w");
    if(itp_ == NULL){
        cout << "Could not open itp file for writing" << endl;
        exit(-1);
    }
    fprintf(itp_, "%s", header_.c_str());
    time_t now = time(0);
    char *dt = ctime(&now);
    fprintf(itp_, "; %s;\n", dt);

    newSection("moleculetype");
    fprintf(itp_, ";molecule name  nrexcl\n");
    fprintf(itp_, "%s %12i\n", resName_.c_str(), 1);
}

ITPWriter::~ITPWriter(){
    if(itp_ != NULL) std::fclose(itp_);
}

void ITPWriter::newSection(const string &section_name){
    section_ = section_name;
    fprintf(itp_, "\n[ %s ]\n", section_.c_str());
}

void ITPWriter::printAtoms(const CGMap &map, const bool isMartini){
    newSection("atoms");
    fprintf(itp_, ";  num  bead type  resnr  resname  bead  cg nr    charge    mass\n");
    double charge = 0.;
    for(BeadMap bead : map.mapping_){
        // MARTINI only has charge on 'Q' beads
        if(isMartini){
            if(bead.type[0] == 'Q'){
                charge = int(bead.charge);
            }else{
                charge = 0.;
            }
        }else{
            charge = bead.charge;
        }

        // Print the line to itp
        fprintf(itp_, "%6i %10s %6i %6s %6s %6i %10.4f",
                bead.num+1, bead.type.c_str(), 1, map.resName_.c_str(),
                bead.name.c_str(), bead.num+1, charge);

        // MARTINI doesn't include masses - all beads are assumed same mass
        if(isMartini){
            fprintf(itp_, ";\n");
        }else{
            fprintf(itp_, " %10.4f;\n", bead.mass);
        }
    }
}

void ITPWriter::printBonds(const BondSet &bond_set){
    newSection("bonds");
    fprintf(itp_, ";atom1 atom2  type equilibrium  force const; unimodality\n");
    for(const BondStruct &bond : bond_set.bonds_){
        // I don't know what the last integer is
        fprintf(itp_, "%5i %5i %5i %12.5f %12.5f; %5.3f\n",
                bond.atomNums_[0]+1, bond.atomNums_[1]+1, 1,
                bond.avg_, bond.forceConstant_, bond.rsqr_);
    }
    newSection("angles");
    fprintf(itp_, ";atom1 atom2 atom3  type equilibrium  force const; unimodality\n");
    for(const BondStruct &bond : bond_set.angles_){
        fprintf(itp_, "%5i %5i %5i %5i %12.5f %12.5f; %5.3f\n",
                bond.atomNums_[0]+1, bond.atomNums_[1]+1,
                bond.atomNums_[2]+1, 2,
                bond.avg_, bond.forceConstant_, bond.rsqr_);
    }
    newSection("dihedrals");
    fprintf(itp_, ";atom1 atom2 atom3 atom4  type equilibrium  force const  mult; unimodality\n");
    for(const BondStruct &bond : bond_set.dihedrals_){
        fprintf(itp_, "%5i %5i %5i %5i %5i %12.5f %12.5f %5i; %5.3f\n",
                bond.atomNums_[0]+1, bond.atomNums_[1]+1,
                bond.atomNums_[2]+1, bond.atomNums_[3]+1,
                1, bond.avg_, bond.forceConstant_, 1, bond.rsqr_);
    }
}