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

/**
* \brief Array of floats with safety features.
*
* Holds a 1, 2 or 3 dimensional array of floats.  In safemode (default) array bounds will be checked
* and Python style negative indices for back counting.  Contains some common array functions which
* can operate on the entire array or on sections.
*/
class ArrayFloat{
protected:
    int dimensions_;
    std::vector<int> size_;
    int elems_;
    float* array_;
    bool fast_;
    bool allocated_;
    int sizex_, sizey_, sizez_;

public:
    int appendedRows_ = 0;
    //ArrayFloat(const std::vector<int> size, const bool fast=false);
    ArrayFloat(const int a, const int b, const int c, const bool fast=false);
    ArrayFloat();
    //TODO Destructor doesn't work??  Fix this
    //~ArrayFloat();
    void init(const int a, const int b=1, const int c=1, const bool fast=false);
    //unsigned long operator ()(int i, int j=0, int k=0) const;
    //unsigned long & operator ()(int i, int j=0, int k=0);
    void append(std::vector<float> vec);
    float& operator()(int x);
    float& operator()(int x, int y);
    float& operator()(int x, int y, int z);
    void zero();
    void linspace(const int a, const int b, const int n, const float min, const float max);
    void linspace(const int a, const int n, const float min, const float max);
    void linspace(const int n, const float min, const float max);
};

/** Linspace (similar to numpy) over a 1d array */
void linspace_1d(array_float_1d array, const float min, const float max, const int steps);

/** Zero a 2d float array */
void zero_float_2d(array_float_2d array, const int a, const int b);
/** Zero a 3d float array */
void zero_float_3d(array_float_3d array, const int a, const int b, const int c);

#endif
