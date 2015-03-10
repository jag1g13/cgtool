#ifndef ARRAYS_H_
#define ARRAYS_H_

#include <vector>


/**
* \brief Array of doubles with safety features.
*
* Holds a 1, 2 or 3 dimensional array of doubles.  In safemode (default) array bounds will be checked
* and Python style negative indices for back counting.  Contains some common array functions which
* can operate on the entire array or on sections.
*/
class Array{
protected:
    /** Dimension of array: allows 1,2,3 */
    int dimensions_;
    /** Size of array in each dimension */
    std::vector<int> size_;
    /** Total number of elements in array */
    int elems_;
    /** Pointer to the actual array */
    double* array_;
    /** Ignore safety features?  Default no */
    bool fast_;
    /** Has array_ been allocated yet? */
    bool allocated_ = false;
    /** Size of the array in the x dimension */
    int sizex_;
    /** Size of the array in the y dimension */
    int sizey_;
    /** Size of the array in the z dimension */
    int sizez_;

public:
    /** How many rows have been appended to the array? */
    int appendedRows_ = 0;

    /** \brief Dimensions_ getter */
    int getDimensions() const {return dimensions_;};
    /** \brief Elems_ getter */
    int getElems() const {return elems_;};

    /** \brief Constructor which allocates the array automatically.
    * The array is zeroed after allocation. */
    Array(const int a, const int b=1, const int c=1, const bool fast=false);

    /** Default constructor which doesn't allocate the array automatically */
    Array();

    /** \brief Destructor.  Just calls Array::free() */
    ~Array();

    /** \brief Initialise the array after calling the default blank constructor.
    * The array is zeroed after allocation. */
    void init(const int a, const int b=1, const int c=1, const bool fast=false);

    /** \brief Append a row to the array into the next blank row.
    * Doesn't have to fill the row. */
    void append(const std::vector<double> &vec);

    /** \brief Append a row to the array into the next blank row.
    * Copies len doubles from *vec.  Doesn't have to fill the row.*/
    void append(const double *vec, const int len);

    /** 1 dimensional access to the array */
    double& operator()(int x);
//    /** 1 dimensional access to the array allowing slicing */
//    double* operator()(int x);
    /** 2 dimensional access to the array */
    double& operator()(int x, int y);
    /** 3 dimensional access to the array */
    double& operator()(int x, int y, int z);

    /** Set all elements to 0.f */
    void zero();

    /** Print all elements of the array */
    void print(const int width=8, const int prec=4, const double scale=1);

    /** Free the array and mark as unallocated */
    void free();

    /** Linspace a line of a 3d array */
    void linspace(const int a, const int b, const int n, const double min, const double max);
    /** Linspace a line of a 2d array */
    void linspace(const int a, const int n, const double min, const double max);
    /** Linspace a line of a 1d array */
    void linspace(const int n, const double min, const double max);

    /** \brief Sum all elements in the array. */
    double sum();

    // operators and friends
    /** Array equality test */
    friend bool operator==(const Array &a, const Array &b);
    /** Array in place subtract operator */
    Array& operator-=(const Array &other);
    /** Array in place add operator */
    Array& operator+=(const Array &other);

    /** RMS difference between two arrays */
    friend double rms(const Array *a, const Array *b);
};

struct StatsBox{
    /** RMS difference between items */
    double rms = 0.f;
    /** RRMS difference between items - RMS / mean value */
    double rrms = 0.f;
    /** Mean difference between items */
    double mean_diff = 0.f;
    /** Mean of values in array A */
    double mean_a = 0.f;
    /** Mean of values in array B */
    double mean_b = 0.f;
};

double vector_rms(const std::vector<double> *a, const std::vector<double> *b);
StatsBox vector_stats(const std::vector<double> *a, const std::vector<double> *b);

#endif
