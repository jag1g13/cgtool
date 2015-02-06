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

    void invert();
    void binHistogram(const int bins);
};

#endif