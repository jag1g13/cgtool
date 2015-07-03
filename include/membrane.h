//
// Created by james on 28/04/15.
//

#ifndef CGTOOL_MEMBRANE_H
#define CGTOOL_MEMBRANE_H

#include <string>
#include <vector>
#include <memory>

#include "frame.h"
#include "array.h"
#include "residue.h"

template <typename T> class LightArray{
protected:
    T *array_;
    int size_[2] = {0, 0};
    int length_ = 0;
public:
    LightArray<T>(){};
    LightArray<T>(const int x, const int y){alloc(x, y);};
    ~LightArray<T>(){
        delete array_;
    };

    /** \brief Copy constructor */
    LightArray<T>(const LightArray &other){
        alloc(other.size_[0], other.size_[1]);
        for(int i=0; i<length_; i++){
            array_[i] = other.array_[i];
        }
    }

    /** \brief Assignment operator */
    LightArray<T>& operator=(const LightArray &other){
       if(size_[0]==other.size_[0] && size_[1]==other.size_[1]){
           for(int i=0; i<length_; i++){
               array_[i] = other.array_[i];
           }
       }
    }

    void alloc(const int x, const int y){
        size_[0] = x; size_[1] = y;
        length_ = x * y;
        array_ = new T[length_];
        if(array_ == nullptr) throw std::runtime_error("Could not allocate array");
    }

    T& operator()(const int x, const int y){
        return array_[x*size_[0] + y];
    }

    const T& at(const int x, const int y) const{
        return array_[x*size_[0] + y];
    }

    /** \brief Apply n iterations of Jacobi smoothing */
    void smooth(const int n_iter=1){
        int jsw = 0;
        int isw = 0;
        for(int ipass = 0; ipass < 2 * n_iter; ipass++){
            jsw = isw;
            for(int i = 1; i < size_[0]-1; i++){
                for(int j = jsw+1; j < size_[1]-1; j += 2){
                    const int loc = i*size_[0] + j;
                    array_[loc] += 0.25 * (array_[loc+size_[0]] + array_[loc+1] + array_[loc-1] + array_[loc-size_[0]] - 4 * array_[loc]);
                }
                jsw = 1 - jsw;
            }
            isw = 1 - isw;
        }
    }

    /** \brief Print the array - format required to account for different types */
    void print(const char *format){
        for(int i=0; i<size_[0]; i++){
            for(int j=0; j<size_[1]; j++){
                printf(format, array_[i*size_[0] + j]);
            }
            printf("\n");
        }
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
    /** Local membrane mean curvature */
    LightArray<double> curvMean_;
    /** Local membrane Gaussian curvature */
    LightArray<double> curvGaussian_;

    /** List of residues present in simulation */
    const std::vector<Residue> *residues_;
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
    FILE *avgFile_ = nullptr;
    /** \brief Number of grid points for each residue */
    std::vector<int> residuePPL_;

    /** \brief Create closest pairs of reference groups between layers */
    void makePairs(const Frame &frame, const std::vector<int> &ref,
                   const std::vector<int> &other, std::map<int, double> &pairs);

    /** \brief Find closest head group to each grid cell */
    double closestLipid(const Frame &frame, const std::vector<int> &ref,
                      const std::map<int, double> &pairs, LightArray<int> &closest);

    void printCSVAreaPerLipid(const float time) const;
    void prepCSVAreaPerLipid();
    void prepCSVAvgThickness();

public:

    /** \brief Print header in CSV or not */
    bool header_ = false;

    /** \brief Blank constructor */
    Membrane(){};

    /** \brief Destructor */
    ~Membrane();

    /** \brief Construct Membrane with vector of Residues present in simulation */
    Membrane(const std::vector<Residue> *residues);

    /** \brief Sort head groups into upper and lower bilayer
     *  Divided into blocks to account for curvature. Size blocks * blocks */
    void sortBilayer(const Frame &frame, const int blocks=4);

    /** \brief Calculate thickness of bilayer */
    void thickness(const Frame &frame, const bool with_reset=false);

    /** \brief Calculate surface area per lipid by residue */
    void areaPerLipid(const LightArray<int> &closest);

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

    /** \brief Set resolution of calculation
     * Number of grid points in x and y */
    void setResolution(const int n);

    /** \brief Zero the running count of membrane thickness */
    void reset();
};

#endif //CGTOOL_MEMBRANE_H
