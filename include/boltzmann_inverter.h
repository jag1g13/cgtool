#ifndef BOLTZMANN_INVERTER_H_
#define BOLTZMANN_INVERTER_H_

#include <vector>

#include "bond_struct.h"
#include "light_array.h"
#include "histogram.h"

/** \brief Class to perform Boltzmann Inversion
*
*/
class BoltzmannInverter{
protected:
    const double temp_;
    const int bins_;

    BondType type_;

    int n_ = 0, meanBin_=0;
    double min_, max_, step_;
    double integral_, mean_, adev_, var_, sdev_;

    /** Store histogram frequencies */
    Histogram histogram_;
    LightArray<double> gaussian_;
    LightArray<double> harmonic_;

    /** \brief Perform a Boltzmann Inversion on a single bond parameter */
    double invertGaussian();
    double invertGaussianSimple();

    /** \brief Sort bond time series into histogram bins */
    void binHistogram(const std::vector<double> &vec);

    /** \brief Calculate R^2 value for calculated gaussian relative to histogram */
    double gaussianRSquared();


public:
    BoltzmannInverter(const double temp=310., const int bins=55);

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