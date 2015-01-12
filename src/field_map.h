#ifndef FIELDMAP_H_
#define FIELDMAP_H_

#include <vector>

#include "arrays.h"
#include "frame.h"
#include "cg_map.h"
//#include "boost/multi_array.hpp"

/**
* \brief Class to hold and operate on electric field maps
*
* Calculates an electric field from point charges in a molecule,
* then fits a set of dipoles to this
*/
class FieldMap{
private:
    /** Dimensions of the field grids (3 ints) */
    std::vector<int> gridDims_;
    /** Array to hold electric field calculated from monopoles */
    ArrayFloat fieldMonopole_;
    /** Array to hold electric field calculated from dipoles */
    ArrayFloat fieldDipole_;
    /** \brief Border to leave around molecule in field grid
    * Also the radius of selection for the CHELPG style grid */
    float border_ = 2.f;    // 2nm
    /** Dipole of each atom, coords, vector, magnitude */
    ArrayFloat dipoles_;
    /** Array to hold grid bounds; needs to be reset each frame (or often) */
    ArrayFloat gridBounds_;
    /** Coordinates of each grid point */
    ArrayFloat gridCoords_;
    /** Centre of grid */
    std::vector<float> gridCentre_;
    ArrayFloat gridContracted_;
    int numGridPoints_;
    std::vector<float> fieldMonopoleContracted_;
    std::vector<float> fieldDipoleContracted_;

public:
    /** Constructor for a blank instance of an electric field map */
    FieldMap();
    /** Constructor for FieldMap to perform setup */
    FieldMap(const int a, const int b, const int c, const int natoms=0);
    /** Determine grid bounds from a Frame object and do setup each time */
    void setupGrid(const Frame *frame);
    /** Calculate the electric field from point charges */
    void calcFieldMonopoles(const Frame *frame);
    /** Calculate the electric field from point dipoles */
    void calcFieldDipoles(const Frame *frame);
    /** \brief Calculate the square of the distance between two points */
    float distSqr(const float *coords, const float x, const float y, const float z);

    /** Create a CHELPG style grid using only points in a shell around the molecule */
    void setupGridContracted(const Frame *frame);
    /** Calculate the electric field from point charges using the CHELPG style grid */
    void calcFieldMonopolesContracted(const Frame *frame);
    /** Calculate the electric field from point dipoles using the CHELPG style grid */
    void calcFieldDipolesContracted(const Frame *frame);
    /** \brief Calculate dipoles on beads directly from frame
    * Modifies charges in atomistic Frame */
    void calcDipolesDirect(const CGMap *cgmap, const Frame *frame, Frame *aa_frame);
    /** Print the dipole array */
    void printDipoles();
    /** Print the electric fields */
    void printFields();
};

/** 3d vector dot product */
inline float dot(const float *A, const float *B);
/** 3d vector magnitude */
inline float abs(const float* vec);

#endif