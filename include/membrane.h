//
// Created by james on 28/04/15.
//

#ifndef CGTOOL_MEMBRANE_H
#define CGTOOL_MEMBRANE_H

#include <string>
#include <vector>

#include "frame.h"
#include "array.h"

struct Res{
    std::string resname;
    std::string ref_atom;
    int num_atoms;
    int num_residues;
};

class Membrane{
protected:
    /** Head group reference atoms in the upper layer */
    std::vector<int> upperHeads_;
    /** Head group reference atoms in the lower layer */
    std::vector<int> lowerHeads_;

    /** The residue that's present - to be made plural soon */
    Res residue_;
    /** Size of the simulation box - assume orthorhombic */
    double box_[3];
    /** 2d Array to hold the membrane thickness on a grid */
    Array thickness_;
    /** 2d Array to hold counts of lipids on a grid */
    Array counts_;
    /** Number of grid points in x and y direction */
    int grid_;
    /** Running average of membrane thickness */
    double avgDist_ = 0.;
    /** Total number of frames processed - for averaging */
    int numFrames_ = 0;

public:

    /** \brief Blank constructor */
    Membrane(){};

    /** \brief Constructor to setup residue */
    Membrane(const std::string &resname, const std::string &ref_atom,
             const int num_atoms, const int num_residues,
             const int grid=20);

    /** \brief Sort head groups into upper and lower bilayer */
    void sortBilayer(const Frame &frame, const int ref_atom);

    /** \brief Calculate thickness of bilayer */
    double thickness(const Frame &frame, const bool with_reset=false);

    /** \brief Calculate thickness with reference to upper or lower leaflet */
    double thickness_with_ref(const Frame &frame, const std::vector<int> &ref,
                              const std::vector<int> &other);

    /** \brief Calculate averages - to be used after all frames have been processed */
    double mean();

    /** \brief Print thickness array to CSV */
    void print_csv(const std::string &filename);
};

#endif //CGTOOL_MEMBRANE_H
