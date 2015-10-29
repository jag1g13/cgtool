//
// Created by james on 28/04/15.
//

#ifndef CGTOOL_MEMBRANE_H
#define CGTOOL_MEMBRANE_H

#include <string>
#include <vector>
#include <set>
#include <array>

#include "frame.h"
#include "residue.h"
#include "light_array.h"

class Membrane{
protected:
    /** Head group reference atoms in the upper layer */
    std::set<int> upperHeads_;
    /** Head group reference atoms in the lower layer */
    std::set<int> lowerHeads_;
    /** Distance from ref in upper leaflet to closest in lower */
    std::map<int, double> upperPair_;
    /** Distance from ref in lower leaflet to closest in upper */
    std::map<int, double> lowerPair_;
    /** Closest lipid to grid point in upper leaflet */
    LightArray<int> closestUpper_;
    /** Closest lipid to grid point in lower leaflet */
    LightArray<int> closestLower_;
    /** Local membrane mean curvature */
    LightArray<double> curvMean_;
    /** Local membrane Gaussian curvature */
    LightArray<double> curvGaussian_;

    /** List of residues present in simulation */
    const std::vector<Residue> &residues_;
    /** Total number of lipids */
    int numLipids_ = 0;
    /** Size of he simulation box - assume orthorhombic */
    std::array<double, 3> box_;
    /** Distance between grid points in xy plane */
    std::array<double, 2> step_;
    /** 2d Array to hold the membrane thickness on a grid */
    LightArray<double> thickness_;
    /** Number of grid points in x and y direction */
    int grid_ = 0;
    /** Total number of frames processed - for averaging */
    int numFrames_ = 0;

    /** \brief File for printing area per lipid */
    FILE *aplFile_ = nullptr;
    FILE *avgFile_ = nullptr;
    /** \brief Number of grid points for each residue */
    std::map<std::string, int> upperResPPL_;
    std::map<std::string, int> lowerResPPL_;
    std::map<std::string, int> upperNumRes_;
    std::map<std::string, int> lowerNumRes_;

    /** \brief Create closest pairs of reference groups between layers */
    void makePairs(const Frame &frame, const std::set<int> &ref,
                   const std::set<int> &other, std::map<int, double> &pairs);

    /** \brief Find closest head group to each grid cell */
    double closestLipid(const Frame &frame, const std::set<int> &ref,
                        const std::map<int, double> &pairs,
                        std::map<std::string, int> &resPPL, LightArray<int> &closest);

    void prepCSVAvgThickness();
    void prepCSVAreaPerLipid();

public:

    /** \brief Print header in CSV or not */
    const bool header_;

    /** \brief Destructor */
    ~Membrane();

    /** \brief Construct Membrane with vector of Residues present in simulation */
    Membrane(const std::vector<Residue> &residues, const Frame &frame,
             const int resolution=100, const int blocks=4, const bool header=true);

    /** \brief Sort head groups into upper and lower bilayer
     *  Divided into blocks to account for curvature. Size blocks * blocks */
    void sortBilayer(const Frame &frame, const int blocks=4);

    /** \brief Calculate thickness of bilayer */
    double thickness(const Frame &frame, const bool with_reset=false);

    /** \brief Calculate curvature of membrane by 2nd order finite differences */
    void curvature(const Frame &frame);

    /** \brief Print curvature grid to file */
    void printCSVCurvature(const std::string &filename) const;

    /** \brief Calculate average thickness */
    double mean() const;

    /** \brief Normalize membrane thickness array in place */
    void normalize(const int smooth_iter=1);

    /** \brief Print thickness array to CSV */
    void printCSV(const std::string &filename) const;

    void printCSVAreaPerLipid(const float time) const;

    /** \brief Set resolution of calculation
     * Number of grid points in x and y */
    void setResolution(const int n);

    /** \brief Zero the running count of membrane thickness */
    void reset();
};

#endif //CGTOOL_MEMBRANE_H
