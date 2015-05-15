//
// Created by james on 07/05/15.
//

#ifndef CGTOOL_RESIDUE_H
#define CGTOOL_RESIDUE_H

#include <string>

#include <assert.h>

class Residue{
public:
    std::string resname = "";
    std::string ref_atom_name = "";
    int ref_atom = -1;
    int num_atoms = -1;
    int num_residues = -1;
    int total_atoms = -1;
    int start = -1;
    int end = -1;
    bool populated = false;

    Residue& operator=(const Residue other);

    void calc_total();

    void init();

    void set_num_atoms(const int val);

    void set_num_residues(const int val);

    void set_start(const int val);

    void set_resname(const std::string &val);

    void set_total_atoms(const int val);

    void print();
};

#endif //CGTOOL_RESIDUE_H
