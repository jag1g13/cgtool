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

struct HeadPair{
    int upper;
    int lower;
};

class Membrane{
protected:
    /** Head group reference atoms in the upper layer */
    std::vector<int> upperHeads_;
    /** Head group reference atoms in the lower layer */
    std::vector<int> lowerHeads_;
    /** Closest in lower leaflet to lipid in upper leaflet */
    std::map<int, int> upperPair_;
    /** Closest in upper leaflet to lipid in lower leaflet */
    std::map<int, int> lowerPair_;

    /** The residue that's present - to be made plural soon */
    Res residue_;
    /** Size of the simulation box - assume orthorhombic */
    double box_[3];
    /** 2d Array to hold the membrane thickness on a grid */
    Array thickness_;
    /** Number of grid points in x and y direction */
    int grid_;
    /** Total number of frames processed - for averaging */
    int numFrames_ = 0;

    /** \brief Create closest pairs of reference groups between layers */
    void makePairs(const Frame &frame, const std::vector<int> &ref,
                   const std::vector<int> &other, std::map<int, int> &pairs);

    /** \brief Calculate thickness with reference to upper or lower leaflet */
    void thicknessWithRef(const Frame &frame, const std::vector<int> &ref,
                            const std::vector<int> &other,
                            const std::map<int, int> &pairs);

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
    void thickness(const Frame &frame, const bool with_reset=false);

    /** \brief Calculate averages - to be used after all frames have been processed */
    double mean();

    /** \brief Normalize membrane thickness array in place */
    void normalize();

    /** \brief Print thickness array to CSV */
    void printCSV(const std::string &filename);
};

#endif //CGTOOL_MEMBRANE_H
