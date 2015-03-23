//
// Created by james on 20/03/15.
//

#ifndef _CGTOOL_SMALL_FUNCTIONS_H_
#define _CGTOOL_SMALL_FUNCTIONS_H_

#include <string>
#include <ctime>

/** \brief Check if a file exists */
bool file_exists(const std::string name);

/** \brief Check if a file exists, if so, rename it.
*
* Makes sure we're not overwriting any existing file.
* Returns true if it's safe to continue */
bool backup_old_file(const std::string name);

/** \brief Print dividers in the text output of a program. */
void split_text_output(const std::string, const clock_t, const int num_threads=1);

/** \brief Dot product of 3d vectors as double[3] */
double dot(const double A[3], const double B[3]);

/** \brief Magnitude of 3d vector as double[3] */
double abs(const double vec[3]);

/** \brief Distance squared between two points as double[3] */
double distSqr(const double c1[3], const double c2[3]);

/** \brief Convert 3d vector as double[3] to polar coordinates */
void polar(const double cart[3], double polar[3]);

#endif //_CGTOOL_SMALL_FUNCTIONS_H_
