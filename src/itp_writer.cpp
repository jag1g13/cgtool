#include "itp_writer.h"

#include <iostream>
#include <cmath>

#include <sysexits.h>

#include "small_functions.h"

using std::string;
using std::fprintf;
using std::cout;
using std::endl;
using std::vector;

ITPWriter::ITPWriter(const string &resName, string itpname, FileFormat format){
    format_ = format;

    char comment = ';';
    switch(format_){
        case FileFormat::GROMACS:
            if(itpname == "") itpname = resName + ".itp";
            comment = ';';
            break;
        case FileFormat::LAMMPS:
            if(itpname == "") itpname = "forcefield." + resName;
            comment = '#';
            break;
    }
    name_ = itpname;
    resName_ = resName;
    backup_old_file(name_.c_str());

    itp_ = std::fopen(name_.c_str(), "w");
    if(itp_ == NULL){
        cout << "Could not open itp file for writing" << endl;
        exit(EX_CANTCREAT);
    }

    // Set the comment marker for this format
    fprintf(itp_, "%c\n", comment);
    for(const string& line : header_) fprintf(itp_, "%c %s", comment, line.c_str());
    fprintf(itp_, "%c\n", comment);

    // Would like timestamp, but it conflicts with testing - can't diff a file with timestamp
//    time_t now = time(0);
//    char *dt = ctime(&now);
//    fprintf(itp_, "; %s;\n", dt);

    newSection("moleculetype");
    fprintf(itp_, ";molecule name  nrexcl\n");
    fprintf(itp_, "%s %12i\n", resName_.c_str(), 1);
}

ITPWriter::~ITPWriter(){
    if(itp_ != NULL) std::fclose(itp_);
}

void ITPWriter::newSection(const string &section_name){
    section_ = section_name;
    switch(format_){
        case FileFormat::GROMACS:
            fprintf(itp_, "\n[ %s ]\n", section_.c_str());
            break;
        case FileFormat::LAMMPS:
            fprintf(itp_, "\n#%s\n", section_.c_str());
            break;
    }
}

void ITPWriter::printAtoms(const CGMap &map, const bool isMartini){
    newSection("atoms");
    fprintf(itp_, ";  num  bead type  resnr  resname  bead  cg nr    charge    mass\n");

    for(BeadMap bead : map.mapping_){
        // MARTINI only has charge on 'Qx' beads
        double charge = bead.charge;
        if(isMartini){
            if(bead.type[0] == 'Q'){
                // Convert to integer with rounding - Mac doesn't have round()
                charge = floor(charge + copysign(0.5, charge));
            }else{
                charge = 0.;
            }
        }

        // Print the line to itp
        fprintf(itp_, "%6i %10s %6i %6s %6s %6i %10.4f",
                bead.num+1, bead.type.c_str(), 1, resName_.c_str(),
                bead.name.c_str(), bead.num+1, charge);

        // MARTINI doesn't include masses - all beads are assumed same mass
        if(!isMartini) fprintf(itp_, " %10.4f", bead.mass);
        fprintf(itp_, ";\n");
    }
}

void ITPWriter::printBonds(const BondSet &bond_set, const bool round){
    const double scale = 3.;
    switch(format_){
        case FileFormat::GROMACS:
            newSection("bonds");
            fprintf(itp_, ";atom1 atom2  type equilibrium  force const; unimodality\n");
            for(const BondStruct &bond : bond_set.bonds_){
                double f_const = bond.forceConstant_;
                if(round) f_const = scale * pow(10, floor(log10(f_const)));
                // I don't know what the last integer is
                fprintf(itp_, "%5i %5i %5i %12.5f %12.5f; %5.3f\n",
                        bond.atomNums_[0]+1, bond.atomNums_[1]+1, 1,
                        bond.avg_, f_const, bond.rsqr_);
            }

            newSection("angles");
            fprintf(itp_, ";atom1 atom2 atom3  type equilibrium  force const; unimodality\n");
            for(const BondStruct &bond : bond_set.angles_){
                double f_const = bond.forceConstant_;
                if(round) f_const = scale * pow(10, floor(log10(f_const)));
                fprintf(itp_, "%5i %5i %5i %5i %12.5f %12.5f; %5.3f\n",
                        bond.atomNums_[0]+1, bond.atomNums_[1]+1,
                        bond.atomNums_[2]+1, 2,
                        bond.avg_, f_const, bond.rsqr_);
            }

            newSection("dihedrals");
            fprintf(itp_, ";atom1 atom2 atom3 atom4  type equilibrium  force const  mult; unimodality\n");
            for(const BondStruct &bond : bond_set.dihedrals_){
                double f_const = bond.forceConstant_;
                if(round) f_const = scale * pow(10, floor(log10(f_const)));
                fprintf(itp_, "%5i %5i %5i %5i %5i %12.5f %12.5f %5i; %5.3f\n",
                        bond.atomNums_[0]+1, bond.atomNums_[1]+1,
                        bond.atomNums_[2]+1, bond.atomNums_[3]+1,
                        1, bond.avg_, f_const, 1, bond.rsqr_);
            }
            break;

        case FileFormat::LAMMPS:
            newSection("bonds");
            fprintf(itp_, "#\t\tbond  eqm  f_const# unimodality\n");
            int i = 1;
            for(const BondStruct &bond : bond_set.bonds_){
                double f_const = bond.forceConstant_;
                if(round) f_const = scale * pow(10, floor(log10(f_const)));
                // I don't know what the last integer is
                fprintf(itp_, "bond_coeff %5i %8.3f %8.3f # %8.3f\n",
                        i, bond.avg_, f_const, bond.rsqr_);
                i++;
            }
            fprintf(itp_, "#\t\tbond  eqm  f_const# unimodality\n");

            newSection("angles");
            fprintf(itp_, "#\t\tbond  eqm  f_const# unimodality\n");
            i = 1;
            for(const BondStruct &bond : bond_set.angles_){
                double f_const = bond.forceConstant_;
                if(round) f_const = scale * pow(10, floor(log10(f_const)));
                // I don't know what the last integer is
                fprintf(itp_, "bond_coeff %5i %8.3f %8.3f # %8.3f\n",
                        i, bond.avg_, f_const, bond.rsqr_);
                i++;
            }
            fprintf(itp_, "#\t\tbond  eqm  f_const# unimodality\n");

            newSection("dihedrals");
            fprintf(itp_, "#\t\tbond  eqm  f_const# unimodality\n");
            i = 1;
            for(const BondStruct &bond : bond_set.dihedrals_){
                double f_const = bond.forceConstant_;
                if(round) f_const = scale * pow(10, floor(log10(f_const)));
                // I don't know what the last integer is
                fprintf(itp_, "bond_coeff %5i %8.3f %8.3f # %8.3f\n",
                        i, bond.avg_, f_const, bond.rsqr_);
                i++;
            }

            break;
    }
}