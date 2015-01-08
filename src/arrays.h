#ifndef ARRAYS_H_
#define ARRAYS_H_

//#include <cstdlib>
#include <vector>

/**
* \brief Array of floats with safety features.
*
* Holds a 1, 2 or 3 dimensional array of floats.  In safemode (default) array bounds will be checked
* and Python style negative indices for back counting.  Contains some common array functions which
* can operate on the entire array or on sections.
*/
class ArrayFloat{
protected:
    /** Dimension of array: allows 1,2,3 */
    int dimensions_;
    /** Size of array in each dimension */
    std::vector<int> size_;
    /** Total number of elements in array */
    int elems_;
    /** Pointer to the actual array */
    float* array_;
    /** Ignore safety features?  Default no */
    bool fast_;
    /** Has array_ been allocated yet? */
    bool allocated_;
    /** Size of the array in the x dimension */
    int sizex_;
    /** Size of the array in the y dimension */
    int sizey_;
    /** Size of the array in the z dimension */
    int sizez_;

public:
    /** How many rows have been appended to the array? */
    int appendedRows_ = 0;
    //ArrayFloat(const std::vector<int> size, const bool fast=false);
    /** \brief Constructor which allocates the array automatically.
    * The array is zeroed after allocation. */
    ArrayFloat(const int a, const int b, const int c, const bool fast=false);
    /** Default constructor which doesn't allocate the array automatically */
    ArrayFloat();
    //TODO Destructor doesn't work??  Fix this
    //~ArrayFloat();
    /** \brief Initialise the array after calling the default blank constructor.
    * The array is zeroed after allocation. */
    void init(const int a, const int b=1, const int c=1, const bool fast=false);
    //unsigned long operator ()(int i, int j=0, int k=0) const;
    //unsigned long & operator ()(int i, int j=0, int k=0);
    /** \brief Append a row to the array into the next blank row.
    * Doesn't have to fill the row. */
    void append(std::vector<float> vec);
    /** \brief Append a row to the array into the next blank row.
    * Copies len floats from *vec.  Doesn't have to fill the row.*/
    void append(const float *vec, const int len);
    /** 1 dimensional access to the array */
    float& operator()(int x);
//    /** 1 dimensional access to the array allowing slicing */
//    float* operator()(int x);
    /** 2 dimensional access to the array */
    float& operator()(int x, int y);
    /** 3 dimensional access to the array */
    float& operator()(int x, int y, int z);
    /** Set all elements to 0.f */
    void zero();
    /** Print all elements of the array */
    void print();
    /** Free the array and mark as unallocated */
    void free();
    /** Linspace a line of a 3d array */
    void linspace(const int a, const int b, const int n, const float min, const float max);
    /** Linspace a line of a 2d array */
    void linspace(const int a, const int n, const float min, const float max);
    /** Linspace a line of a 1d array */
    void linspace(const int n, const float min, const float max);

    // operators
    friend bool operator==(ArrayFloat &a, ArrayFloat &b);
};

#endif
