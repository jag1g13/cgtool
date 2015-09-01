//
// Created by james on 14/05/15.
//

#include "residue.h"

#include <cstdio>

#include <assert.h>
#include <sysexits.h>

using std::string;

Residue& Residue::operator=(const Residue other){
    resname = other.resname;
    ref_atom = other.ref_atom;
    num_atoms = other.num_atoms;
    num_residues = other.num_residues;
    total_atoms = other.total_atoms;
    start = other.start;
    end = other.end;

    return *this;
}

void Residue::calc_total(){
    total_atoms = num_atoms * num_residues;
    end = start + total_atoms;
}

void Residue::init(){
    num_atoms = 0;
    num_residues = 0;
}

void Residue::set_num_atoms(const int val){
    if(num_atoms != -1) assert(val == num_atoms);
    num_atoms = val;
}

void Residue::set_num_residues(const int val){
    if(num_residues != -1) assert(val == num_residues);
    num_residues = val;
}

void Residue::set_start(const int val){
//    if(start != -1) assert(val == start);
    start = val;
}

void Residue::set_resname(const string &val){
    if(resname.compare("") && val.compare(resname)){
        // If already set
        printf("Residue %5s in GRO did not match %5s in CFG\n",
               val.c_str(), resname.c_str());
        exit(EX_USAGE);
    }else{
        resname = val;
    }
}

void Residue::set_total_atoms(const int val){
    // Assertion not true if both num_residues and num_atoms are -1
//    if(total_atoms != -1) assert(val == total_atoms);
    total_atoms = val;
}

void Residue::print(const bool extra) const{
    if(extra){
        printf("%'6d x %5s with %'3d atoms total %'6d, starting at %'6d end %'6d\n",
               num_residues, resname.c_str(), num_atoms, total_atoms, start, end);
    }else{

        printf("%'6d x %5s with %'3d atoms\n",
               num_residues, resname.c_str(), num_atoms);
    }
}
