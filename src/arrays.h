#ifndef ARRAYS_H_
#define ARRAYS_H_

//#include <cstdlib>
#include <vector>

typedef float* array_float_1d;
typedef float** array_float_2d;
typedef float*** array_float_3d;

/** Allocate a 2d array of floats, return pointer if successful */
array_float_2d alloc_float_2d(const int a, const int b);
/** Allocate a 3d array of floats, return pointer if successful */
array_float_3d alloc_float_3d(const int a, const int b, const int c);

/** Allocate a 2d array of floats in flat form, return pointer if successful */
array_float_1d alloc_float_2d_flat(const int a, const int b, const int c);
/** Allocate a 3d array of floats in flat form, return pointer if successful */
array_float_1d alloc_float_3d_flat(const int a, const int b, const int c);

class Array{
protected:
    int dimensions_;
    std::vector<int> size_;

public:
    Array(std::vector<int> size);
};

class ArrayFloatOne : Array{

};

/** Linspace (similar to numpy) over a 1d array */
void linspace_1d(array_float_1d array, const float min, const float max, const int steps);

/** Zero a 2d float array */
void zero_float_2d(array_float_2d array, const int a, const int b);
/** Zero a 3d float array */
void zero_float_3d(array_float_3d array, const int a, const int b, const int c);

#endif
