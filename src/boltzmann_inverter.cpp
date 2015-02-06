#include "boltzmann_inverter.h"

#include <iostream>

using std::cout;
using std::endl;

void BoltzmannInverter::invert(){
}

void BoltzmannInverter::binHistogram(const int bins){
    BondStruct *bond;
    float max = float(bond->avg_), min = float(bond->avg_);
    histogram_.init(bins);

    for(const float val : bond->values_){
        if(val < min) min = val;
        if(val > max) max = val;
    }

    float step = (max - min) / (bins-1);

    for(const float val : bond->values_){
        int loc = int((val - min) / step);
        if(loc < 0 || loc > bins-1) cout << loc << endl;
        histogram_(loc)++;
    }
}
