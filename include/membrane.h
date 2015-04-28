//
// Created by james on 28/04/15.
//

#ifndef CGTOOL_MEMBRANE_H
#define CGTOOL_MEMBRANE_H

#include <string>
#include <vector>

#include "frame.h"

struct Res{
    std::string resname;
    std::string ref_atom;
    int num_atoms;
    int num_residues;
};

class Membrane{
protected:
    std::vector<int> upperHeads_;
    std::vector<int> lowerHeads_;
    Res residue_;

public:

    /** \brief Blank constructor */
    Membrane(){};

    /** \brief Constructor to setup residue */
    Membrane(const std::string &resname, const std::string &ref_atom, const int num_atoms, const int num_residues);

    /** \brief Sort head groups into upper and lower bilayer */
    void sortBilayer(const Frame &frame, const int ref_atom);

    /** \brief Calculate thickness of bilayer */
    double thickness(const Frame &frame);
};

#endif //CGTOOL_MEMBRANE_H
