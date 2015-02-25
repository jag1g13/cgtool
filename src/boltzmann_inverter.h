#ifndef BOLTZMANN_INVERTER_H_
#define BOLTZMANN_INVERTER_H_

#include "bondset.h"
#include "array.h"

/** \brief Class to perform Boltzmann Inversion
*
*/
class BoltzmannInverter{
protected:
    /** Pointer to associated BondSet */
    BondSet *bondSet_;
    /** Store histogram frequencies */
    ArrayFloat histogram_;
    ArrayFloat gaussian_;
    int bins_ = 0, n_ = 0;
    double min_, max_, step_;
    double mean_, adev_, var_, sdev_;

public:
    BoltzmannInverter(){};
    /** \brief Construct a new Boltzmann Inverter and assign a bondset */
    BoltzmannInverter(BondSet *bondSet){bondSet_ = bondSet;};

    /** \brief Perform a Boltzmann Inversion on a single bond parameter */
    void invertGaussian();
    /** \brief Sort bond time series into histogram bins */
    void binHistogram(const BondStruct &bond, const int bins=100);

    /** \brief Calculate statistical moments of bond data.
    * Mean, standard deviation, skewness and kurtosis.
    * This data may not be useful for a multi-modal distribution.
    */
    void statisticalMoments(const std::vector<double> &vec);
//    void statisticalMoments(ArrayFloat &array);

    /** \brief Calculate R^2 value for calculated gaussian relative to histogram */
    double gaussianRSquared();
};

#endif