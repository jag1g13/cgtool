#ifndef ARRAYS_H_
#define ARRAYS_H_

#include <vector>
#include <string>


//TODO try out Boost/FLENS arrays, are they better?
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
    int size_[3];
    /** Total number of elements in array */
    int elems_;
    /** Pointer to the actual array */
    double* array_;
    /** Ignore safety features?  Default no */
    bool fast_;
    /** Has array_ been allocated yet? */
    bool allocated_ = false;
    /** How many rows have been appended to the array? */
    int appendedRows_ = 0;

public:
    /** \brief Dimensions_ getter */
    int getDimensions() const {return dimensions_;};
    /** \brief Elems_ getter */
    int getElems() const {return elems_;};
    /** \brief AppendedRows_ getter */
    int getAppendedRows() const {return appendedRows_;};
    /** \brief Set appendedRows_ to zero.  Start appending from beginning. */
    void resetAppendedRows(){appendedRows_ = 0;};

    /** \brief Constructor which allocates the array automatically.
    * The array is zeroed after allocation.  Calls init(). */
    Array(const int a, const int b=1, const int c=1, const bool fast=false);

    /** \brief Empty constructor for Array, does nothing */
    Array(){};

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
    void append(const double *vec, int len=-1);

    /** 1 dimensional access to the array */
    double& operator()(int x);
    /** 2 dimensional access to the array */
    double& operator()(int x, int y);
    /** 3 dimensional access to the array */
    double& operator()(int x, int y, int z);

    /** 2 dimensional access to the array - extra bounds checking */
    double& at(int x, int y);

    /** Set all elements to 0. */
    void zero();

    /** Print all elements of the array */
    void print(const int width=8, const int prec=4, const double scale=1);
    /** Print array to CSV */
    void print_csv(const std::string &filename, const int remove_border=0);

    /** Free the array and mark as unallocated */
    void free();

    /** Linspace a line of a 3d array */
    void linspace(const int a, const int b, const int n, const double min, const double max);
    /** Linspace a line of a 2d array */
    void linspace(const int a, const int n, const double min, const double max);
    /** Linspace a line of a 1d array */
    void linspace(const int n, const double min, const double max);

    /** \brief Sum all elements in the array */
    double sum();
    /** \brief Mean of all elements in array */
    double mean();

    /** \brief Apply smoothing to array
     * Uses n_iter Gauss-Seidel Red-Black iterations */
    void smooth(const int n_iter=1);
    /** \brief Find all zero entries and interpolate values from neighbours */
    void interpolate_zeros();

    // operators and friends
    /** Array equality test */
    friend bool operator==(const Array &a, const Array &b);
    /** Array in place subtract operator */
    Array& operator-=(const Array &other);
    /** Array in place add operator */
    Array& operator+=(const Array &other);

    /** Elementwise in place multiply */
    void element_multiply(const Array &other);
    /** Elementwise in place divide */
    void element_divide(const Array &other);

    /** Replace all NaN entries with the mean value */
    void replace_nan();

    /** RMS difference between two arrays */
    friend double rmsd(const Array &a, const Array &b);
};


#endif
