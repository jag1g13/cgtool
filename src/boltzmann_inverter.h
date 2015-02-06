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

public:
    BoltzmannInverter(){};
    BoltzmannInverter(BondSet *bondSet){bondSet_ = bondSet;};

    /** \brief Perform a Boltzmann Inversion on a single bond parameter */
    void invert();
    /** \brief Sort bond time series into histogram bins */
    void binHistogram(const BondStruct &bond, const int bins);
};

#endif