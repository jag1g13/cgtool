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

ITPWriter::ITPWriter(const vector<Residue> &residues, const FileFormat file_format,
                     const FieldFormat field_format, string itpname){
    format_ = file_format;
    fieldFormat_ = field_format;
    resName_ = residues[0].resname;

    switch(format_){
        case FileFormat::GROMACS:
            if(itpname == "") itpname = resName_ + ".itp";
            comment_ = ';';
            break;

        case FileFormat::LAMMPS:
            // Considering just leaving filenames the same
            if(itpname == "") itpname = "forcefield." + resName_;
            comment_ = '#';
            break;
    }
    name_ = itpname;
    backup_old_file(name_.c_str());

    itp_ = std::fopen(name_.c_str(), "w");
    if(itp_ == NULL){
        cout << "Could not open itp file for writing" << endl;
        exit(EX_CANTCREAT);
    }

    // Set the comment marker for this format
    fprintf(itp_, "%c\n", comment_);
    for(const string& line : header_) fprintf(itp_, "%c %s", comment_, line.c_str());
    fprintf(itp_, "%c\n", comment_);

    // Would like timestamp, but it conflicts with testing - can't diff a file with timestamp
//    time_t now = time(0);
//    char *dt = ctime(&now);
//    fprintf(itp_, "; %s;\n", dt);

    switch(format_){
        case FileFormat::GROMACS:
            newSection("moleculetype");
            fprintf(itp_, ";molecule name  nrexcl\n");
            fprintf(itp_, "%s %12i\n", resName_.c_str(), 1);
            break;
        case FileFormat::LAMMPS:
            fprintf(itp_, "\n#molecule name %s\n", resName_.c_str());
            break;
    }
}

ITPWriter::~ITPWriter(){
    if(itp_ != NULL) std::fclose(itp_);
}

void ITPWriter::newSection(const string &section_name) const{
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

void ITPWriter::printAtoms(const CGMap &map) const{
    switch(format_){
        case FileFormat::GROMACS:
            newSection("atoms");
            fprintf(itp_, ";  num   beadtype  resnr   resnm  bead  chrg#     charge    mass\n");

            for(BeadMap bead : map.mapping_){
                // MARTINI only has charge on 'Qx' beads
                double charge = bead.charge;
                if(fieldFormat_ == FieldFormat::MARTINI){
                    if(bead.type[0] == 'Q'){
                        // Convert to integer with rounding - Mac doesn't have round()
                        charge = floor(charge + copysign(0.5, charge));
                    }else{
                        charge = 0.;
                    }
                }

                fprintf(itp_, "%6i %10s %6i %6s %6s %6i %10.4f",
                        bead.num+1, bead.type.c_str(), 1, resName_.c_str(),
                        bead.name.c_str(), bead.num+1, charge);

                // MARTINI doesn't include masses - all beads are assumed same mass
                if(fieldFormat_ != FieldFormat::MARTINI) fprintf(itp_, " %10.4f", bead.mass);
                fprintf(itp_, ";\n");
            }
            break;

        case FileFormat::LAMMPS:
            fprintf(itp_, "\n#atom masses\n");

            for(BeadMap bead : map.mapping_){
                fprintf(itp_, "mass %4i %8.3f\n",
                        bead.num+1, bead.mass);
            }
            break;
    }
}

void ITPWriter::printBonds(const BondSet &bond_set, const bool round) const{
    const double scale = 3.;
    switch(format_){
        case FileFormat::GROMACS:
            newSection("bonds");
            fprintf(itp_, ";atm1  atm2  type  equilibrium  force const  unimodality\n");
            for(const BondStruct &bond : bond_set.bonds_){
                double f_const = bond.forceConstant_;
                if(round) f_const = scale * pow(10, floor(log10(f_const)));
                fprintf(itp_, "%5i %5i %5i %12.5f %12.5f; %5.3f\n",
                        bond.atomNums_[0]+1, bond.atomNums_[1]+1, 1,
                        bond.avg_, f_const, bond.rsqr_);
            }

            newSection("angles");
            fprintf(itp_, ";atm1  atm2  atm3  type  equilibrium  force const  unimodality\n");
            for(const BondStruct &bond : bond_set.angles_){
                double f_const = bond.forceConstant_;
                if(round) f_const = scale * pow(10, floor(log10(f_const)));
                fprintf(itp_, "%5i %5i %5i %5i %12.5f %12.5f; %5.3f\n",
                        bond.atomNums_[0]+1, bond.atomNums_[1]+1,
                        bond.atomNums_[2]+1, 2,
                        bond.avg_, f_const, bond.rsqr_);
            }

            newSection("dihedrals");
            fprintf(itp_, ";atm1  atm2  atm3  atm4  type  equilibrium  force const  mult  unimodality\n");
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
            fprintf(itp_, "\n#bonds     bond    equil  f_const  unimodality\n");
            int i = 1;
            for(const BondStruct &bond : bond_set.bonds_){
                fprintf(itp_, "bond_coeff %4i %8.3f %8.3f  #  %8.3f\n",
                        i, bond.avg_, bond.forceConstant_, bond.rsqr_);
                i++;
            }

            fprintf(itp_, "\n#angles     bond                    equil  f_const  unimodality\n");
            i = 1;
            for(const BondStruct &bond : bond_set.angles_){
                fprintf(itp_, "angle_coeff %4i  cosine/squared %8.3f %8.3f  #  %8.3f\n",
                        i, bond.avg_, bond.forceConstant_, bond.rsqr_);
                i++;
            }

//            printBondSection(bond_set.dihedrals_);

            break;
    }
}

void ITPWriter::printAtomTypes(const CGMap &cgmap) const{
    switch(format_){
        case FileFormat::GROMACS:
            newSection("atomtype");
            fprintf(itp_, ";name  at.num   mass      charge  ptype       c6           c12\n");
            for(const BeadMap &bead : cgmap.mapping_){
                fprintf(itp_, "%5s%5i%11.3f%11.3f%6s%14.10f%14.10f\n",
                        bead.type.c_str(), 0, bead.mass, bead.charge,
                        "A", bead.c06, bead.c12);
            }
            break;
        case FileFormat::LAMMPS:
            printf("LAMMPS Non-bonded output not yet supported\n");
            exit(EX_USAGE);
    }
}