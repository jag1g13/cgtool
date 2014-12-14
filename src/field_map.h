#ifndef FIELDMAP_H_
#define FIELDMAP_H_

#include <vector>

#include "general.h"
#include "frame.h"
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
    array_float_3d fieldMonopole_;
    /** Array to hold electric field calculated from dipoles */
    array_float_3d fieldDipole_;
    /** Border to leave around molecule in field grid */
    float border_ = 1.f;
    /** Array to hold atomic dipoles */
    array_float_2d dipoles_;
    /** Array to hold grid bounds; needs to be reset each frame (or often) */
    array_float_2d gridBounds_;
    /** Coordinates of each grid point */
    array_float_2d gridCoords_;

public:
    /** Constructor for a blank instance of an electric field map */
    FieldMap();
    /** Constructor for FieldMap to perform setup */
    FieldMap(const int a, const int b, const int c);
    /** Determine grid bounds from a Frame object and do setup each time */
    void setupGrid(Frame *frame);
};

#endif