//
// Created by james on 20/03/15.
//

#ifndef _CGTOOL_SMALL_FUNCTIONS_H_
#define _CGTOOL_SMALL_FUNCTIONS_H_

#include <string>
#include <ctime>
#include <vector>
#include <cmath>

#include <time.h>

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
inline double dot(const double A[3], const double B[3]){
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

/** \brief Magnitude of 3d vector as double[3] */
inline double abs(const double vec[3]){
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

/** \brief Distance squared between two points as double[3] */
inline double distSqr(const double c1[3], const double c2[3]){
    return (c1[0] - c2[0]) * (c1[0] - c2[0]) +
           (c1[1] - c2[1]) * (c1[1] - c2[1]) +
           (c1[2] - c2[2]) * (c1[2] - c2[2]);
}

/** \brief Distance squared between two points in a plane as double[3] */
inline double distSqrPlane(const double c1[3], const double c2[3]){
    return (c1[0] - c2[0]) * (c1[0] - c2[0]) +
           (c1[1] - c2[1]) * (c1[1] - c2[1]);
}

/** \brief Convert 3d vector as double[3] to polar coordinates */
void polar(const double cart[3], double polar[3]);

/** \brief Get (approximate) number of frames in XTC file */
int get_xtc_num_frames(const std::string &xtcname);

struct StatsBox{
    /** RMS difference between items */
    double rmsd = 0.;
    /** RRMS difference between items - RMS / mean value */
    double nrmsd = 0.;
    /** Mean difference between items */
    double diff_means = 0.;
    /** Mean of values in array A */
    double mean_a = 0.;
    /** Mean of values in array B */
    double mean_b = 0.;
    /** Minimum value of a **/
    double min_a;
    /** Maximum value of a **/
    double max_a;
};

StatsBox vector_stats(const std::vector<double> &a, const std::vector<double> &b);

#endif //_CGTOOL_SMALL_FUNCTIONS_H_
