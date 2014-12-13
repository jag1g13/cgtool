#ifndef FIELDMAP_H_
#define FIELDMAP_H_

#include <vector>

#include "boost/multi_array.hpp"

/**
* \brief Class to hold and operate on electric field maps
*
* Calculates an electric field from point charges in a molecule,
* then fits a set of dipoles to this
*/
class FieldMap{
private:
    typedef boost::multi_array<float, 2> array_float_2d;
    typedef array_float_2d::index index;
    /** Dimensions of the field grids (3 ints) */
    std::vector<int> gridDims_;
    /** Array to hold electric field calculated from monopoles */
    array_float_2d fieldMonopole_();
    //std::vector<std::vector<float>> fieldMonopole;
    /** Array to hold electric field calculated from dipoles */
    array_float_2d fieldDipole_;
    //std::vector<std::vector<float>> fieldDipole;
    //std::vector<std::vector<float>> dipoles;
    float border_ = 1.f;
};

#endif