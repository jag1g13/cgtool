//
// Created by james on 28/04/15.
//

#ifndef CGTOOL_MEMBRANE_H
#define CGTOOL_MEMBRANE_H

#include <string>
#include <vector>

#include "frame.h"
#include "array.h"
#include "residue.h"

class Membrane{
protected:
    /** Head group reference atoms in the upper layer */
    std::vector<int> upperHeads_;
    /** Head group reference atoms in the lower layer */
    std::vector<int> lowerHeads_;
    /** Distance from ref in upper leaflet to closest in lower */
    std::map<int, double> upperPair_;
    /** Distance from ref in lower leaflet to closest in upper */
    std::map<int, double> lowerPair_;

    /** The residue that's present - to be made plural soon */
    Residue residue_;
    std::vector<Residue> residues_;
    /** Size of he simulation box - assume orthorhombic */
    double box_[3];
    /** Distance between grid points in xy plane */
    double step_[2];
    /** 2d Array to hold the membrane thickness on a grid */
    Array thickness_;
    /** Number of grid points in x and y direction */
    int grid_;
    /** Total number of frames processed - for averaging */
    int numFrames_ = 0;

    /** \brief Create closest pairs of reference groups between layers */
    void makePairs(const Frame &frame, const std::vector<int> &ref,
                   const std::vector<int> &other, std::map<int, double> &pairs);

    /** \brief Calculate thickness with reference to upper or lower leaflet */
    void thicknessWithRef(const Frame &frame, const std::vector<int> &ref,
                            const std::vector<int> &other,
                            const std::map<int, double> &pairs);

public:

    /** \brief Blank constructor */
    Membrane(){};

    Membrane(const Residue &residue);
    Membrane(const std::vector<Residue> &residues);

    /** \brief Sort head groups into upper and lower bilayer */
    void sortBilayer(const Frame &frame, const int ref_atom);

    /** \brief Calculate thickness of bilayer */
    void thickness(const Frame &frame, const bool with_reset=false);

    /** \brief Calculate averages - to be used after all frames have been processed */
    double mean();

    /** \brief Normalize membrane thickness array in place */
    void normalize(const int smooth_iter=1);

    /** \brief Print thickness array to CSV */
    void printCSV(const std::string &filename);

    /** \brief Set resolution of calculation
     * Number of grid points in x and y */
    void setResolution(const int n);
};

#endif //CGTOOL_MEMBRANE_H
