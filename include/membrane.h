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

template <class type> class LightArray{
protected:
    type *array_ = nullptr;
    int size_[2] = {0, 0};
public:
    LightArray(){};
    LightArray(const int x, const int y){alloc(x, y);};
    ~LightArray(){if(array_ != nullptr) delete[] array_;};

    void alloc(const int x, const int y){
        size_[0] = x; size_[1] = y;
        array_ = new type[x*y];
    }

    type& operator()(const int x, const int y){
        return array_[x*size_[0] + y];
    }

    const type& at(const int x, const int y) const{
        return array_[x*size_[0] + y];
    }


};

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
    /** Closest lipid to grid point in upper leaflet */
    LightArray<int> closestUpper_;
    /** Closest lipid to grid point in lower leaflet */
    LightArray<int> closestLower_;

    /** List of residues present in simulation */
    std::vector<Residue> residues_;
    /** Total number of lipids */
    int numLipids_ = 0;
    /** Size of he simulation box - assume orthorhombic */
    double box_[3];
    /** Distance between grid points in xy plane */
    double step_[2];
    /** 2d Array to hold the membrane thickness on a grid */
    Array thickness_;
    /** Number of grid points in x and y direction */
    int grid_ = 0;
    /** Total number of frames processed - for averaging */
    int numFrames_ = 0;

    /** \brief File for printing area per lipid */
    FILE *aplFile_ = nullptr;
    /** \brief Number of grid points for each residue */
    std::vector<int> residuePPL_;

    /** \brief Create closest pairs of reference groups between layers */
    void makePairs(const Frame &frame, const std::vector<int> &ref,
                   const std::vector<int> &other, std::map<int, double> &pairs);

    /** \brief Find closest head group to each grid cell */
    void closestLipid(const Frame &frame, const std::vector<int> &ref,
                      const std::map<int, double> &pairs, LightArray<int> &closest);

    void printCSVAreaPerLipid(const float time) const;
    void prepCSVAreaPerLipid();

public:

    /** \brief Print header in CSV or not */
    bool header_ = false;

    /** \brief Blank constructor */
    Membrane(){};

    /** \brief Destructor */
    ~Membrane(){fclose(aplFile_);};

    /** \brief Construct Membrane with vector of Residues present in simulation */
    Membrane(const std::vector<Residue> &residues);

    /** \brief Sort head groups into upper and lower bilayer
     *  Divided into blocks to account for curvature. Size blocks * blocks */
    void sortBilayer(const Frame &frame, const int blocks=4);

    /** \brief Calculate thickness of bilayer */
    void thickness(const Frame &frame, const bool with_reset=false);

    /** \brief Calculate surface area per lipid by residue */
    void areaPerLipid(const LightArray<int> &closest);

    /** \brief Calculate curvature of membrane by 2nd order finite differences */
    void curvature(const LightArray<int> &upper, const LightArray<int> &lower,
                   const Frame &frame);

    /** \brief Calculate average thickness */
    double mean() const;

    /** \brief Normalize membrane thickness array in place */
    void normalize(const int smooth_iter=1);

    /** \brief Print thickness array to CSV */
    void printCSV(const std::string &filename) const;

    /** \brief Set resolution of calculation
     * Number of grid points in x and y */
    void setResolution(const int n);

    /** \brief Zero the running count of membrane thickness */
    void reset();
};

#endif //CGTOOL_MEMBRANE_H
