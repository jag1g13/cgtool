#ifndef BOLTZMANN_INVERTER_H_
#define BOLTZMANN_INVERTER_H_

#include <vector>

#include "bond_struct.h"
#include "array.h"
#include "histogram.h"

/** \brief Class to perform Boltzmann Inversion
*
*/
class BoltzmannInverter{
protected:
    /** Store histogram frequencies */
    Histogram histogram_;
    Array gaussian_;
    Array harmonic_;
    int bins_ = 55, n_ = 0, meanBin_=0;
    double temp_ = 310.;
    double min_, max_, step_;
    double integral_, mean_, adev_, var_, sdev_;
    BondType type_ = BondType::LENGTH;

    /** \brief Print an array/histogram to terminal for debugging */
    void printGraph(Array &arr, const int scale=10);

    /** \brief Perform a Boltzmann Inversion on a single bond parameter */
    double invertGaussian();
    double invertGaussianSimple();

    /** \brief Sort bond time series into histogram bins */
    void binHistogram(const std::vector<double> &vec);

    /** \brief Calculate R^2 value for calculated gaussian relative to histogram */
    double gaussianRSquared();


public:
    BoltzmannInverter(const double temp=310, const int bins=-1);

    /** \brief Calculate statistical moments of bond data.
    * Mean, standard deviation, skewness and kurtosis.
    * This data may not be useful for a multi-modal distribution.
     * \returns Mean of vector
    */
    double statisticalMoments(const std::vector<double> &vec);

    /** Perform all of the necessary calculations to get a force constant */
    void calculate(BondStruct &bond);
};

#endif