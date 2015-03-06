#ifndef BOLTZMANN_INVERTER_H_
#define BOLTZMANN_INVERTER_H_

#include "bondset.h"
#include "array.h"

typedef unsigned int uint;

/** \brief Class to perform Boltzmann Inversion
*
*/
class BoltzmannInverter{
protected:
    /** Store histogram frequencies */
    Array histogram_;
    Array gaussian_;
    Array harmonic_;
    int bins_ = 0, n_ = 0;
    double min_, max_, step_, maxGau_=0.;
    double integral_, mean_, adev_, var_, sdev_;

    /** \brief Print an array/histogram to terminal for debugging */
    void printGraph(Array &arr, const int scale=10);

public:
    BoltzmannInverter(){};

    /** \brief Perform a Boltzmann Inversion on a single bond parameter */
    double invertGaussian();
    /** \brief Sort bond time series into histogram bins */
    void binHistogram(const BondStruct &bond, const int bins=100);

    /** \brief Calculate statistical moments of bond data.
    * Mean, standard deviation, skewness and kurtosis.
    * This data may not be useful for a multi-modal distribution.
    */
    void statisticalMoments(const std::vector<double> &vec);

    /** \brief Calculate R^2 value for calculated gaussian relative to histogram */
    double gaussianRSquared();
};

#endif