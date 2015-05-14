//
// Created by james on 07/05/15.
//

#ifndef CGTOOL_RESIDUE_H
#define CGTOOL_RESIDUE_H

#include <assert.h>

struct Residue{
    std::string resname = "";
    std::string ref_atom_name = "";
    int ref_atom = -1;
    int num_atoms = -1;
    int num_residues = -1;
    int total_atoms = -1;
    int start = -1;
    bool populated = false;
    Residue& operator=(const Residue other);

    void calc_total(){
        total_atoms = num_atoms * num_residues;
    };

    void init(){
        num_atoms = 0;
        num_residues = 0;
    }

    void set_num_atoms(const int val){
        if(num_atoms != -1) assert(val == num_atoms);
        num_atoms = val;
    };

    void set_num_residues(const int val){
        if(num_residues != -1) assert(val == num_residues);
        num_residues = val;
    };

    void set_start(const int val){
        if(start != -1) assert(val == start);
        start = val;
    };

    void set_resname(const std::string &val){
        if(resname.compare("")) assert(!val.compare(resname));
        resname = val;
    };

    void set_total_atoms(const int val){
        if(total_atoms != -1) assert(val == total_atoms);
        total_atoms = val;
    };

    void print(){
        printf("%'6d x %5s with %'3d atoms, starting at %'6d\n",
               num_residues, resname.c_str(), num_atoms, start);
    };
};

#endif //CGTOOL_RESIDUE_H
