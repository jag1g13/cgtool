//
// Created by james on 20/03/15.
//

#ifndef _CGTOOL_SMALL_FUNCTIONS_H_
#define _CGTOOL_SMALL_FUNCTIONS_H_

#include <string>
#include <ctime>
#include <vector>
#include <cmath>
#include <array>
#include <cassert>


/** \brief Check if a file exists */
bool file_exists(const std::string name);

/** \brief Get size of file */
long file_size(const std::string filename);

/** \brief Get clock time in seconds - thread stable */
double start_timer();

/** \brief How many seconds have passed since the given time? */
double end_timer(const double since);

/** \brief Check if a file exists, if so, rename it.
*
* Makes sure we're not overwriting any existing file.
* Returns true if it's safe to continue */
bool backup_old_file(const std::string name);

/** \brief Print dividers in the text output of a program. */
void split_text_output(const std::string &name, const double start);

/** \brief Dot product of 3d vectors as double[3] */
template<std::size_t SIZE>
inline double dot(const std::array<double, SIZE> &A, const std::array<double, SIZE> &B){
    assert(SIZE >= 3);
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

template<std::size_t SIZE>
inline void cross(const std::array<double, SIZE> &A, const std::array<double, SIZE> &B,
                  std::array<double, SIZE> &C){
    assert(SIZE >= 3);
    C[0] = A[1]*B[2] - A[2]*B[1];
    C[1] = A[0]*B[2] - A[2]*B[0];
    C[2] = A[0]*B[1] - A[1]*B[0];
}

/** \brief Magnitude of 3d vector as std::array<double, 3>*/
template<std::size_t SIZE>
inline double abs(const std::array<double, SIZE> &vec){
    assert(SIZE >= 3);
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

template<std::size_t SIZE>
std::array<double, SIZE> operator-(const std::array<double, SIZE> &vec,
                                   const std::array<double, SIZE> &vec2){
    std::array<double, SIZE> res;
    for(std::size_t i=0; i<SIZE; i++) res[i] = vec[i] - vec2[i];
    return res;
}

/** \brief Distance squared between two points as std::array<double, 3> */
template<std::size_t SIZE>
inline double distSqr(const std::array<double, SIZE> &c1, const std::array<double, SIZE> &c2){
    assert(SIZE >= 3);
    return (c1[0] - c2[0]) * (c1[0] - c2[0]) +
           (c1[1] - c2[1]) * (c1[1] - c2[1]) +
           (c1[2] - c2[2]) * (c1[2] - c2[2]);
}

/** \brief Distance squared between two points in a plane as std::array<double, 3>
 *   Slightly more efficient than distSqr */
template<std::size_t SIZE>
inline double distSqrPlane(const std::array<double, SIZE> &c1, const std::array<double, SIZE> &c2){
    assert(SIZE >= 2);
    return (c1[0] - c2[0]) * (c1[0] - c2[0]) +
           (c1[1] - c2[1]) * (c1[1] - c2[1]);
}


/** \brief Convert 3d vector as double[3] to polar coordinates */
void polar(const double cart[3], double polar[3]);

/** \brief Get number of frames in XTC file */
int get_xtc_num_frames(const std::string &xtcname);

/** \brief Calculate mean of vector */
double vector_mean(std::vector<double> &vec);

/** \brief Calculate standard error of vector using Welford's one pass method */
double vector_stderr(std::vector<double> &vec);

#endif //_CGTOOL_SMALL_FUNCTIONS_H_
