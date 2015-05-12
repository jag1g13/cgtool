//
// Created by james on 07/05/15.
//

#ifndef CGTOOL_RESIDUE_H
#define CGTOOL_RESIDUE_H

struct Residue{
    std::string resname = "";
    std::string ref_atom_name = "";
    int ref_atom = 0;
    int num_atoms = -1;
    int num_residues = -1;
    int total_atoms = -1;
    int start = 0;
    Residue& operator=(const Residue other);
    void calc_total(){total_atoms = num_atoms * num_residues;};
};

#endif //CGTOOL_RESIDUE_H
