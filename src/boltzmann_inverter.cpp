#include "boltzmann_inverter.h"

#include <iostream>

#include <math.h>

using std::cout;
using std::endl;

void BoltzmannInverter::invertGaussian(){
    /* line from Python version
    y_inv = -R * T * np.log(y_fit / (x_fit*x_fit))
     */
}

void BoltzmannInverter::binHistogram(const BondStruct &bond, const int bins){
    double max = bond.avg_, min = bond.avg_;
    histogram_.init(bins);

    for(const double val : bond.values_){
        if(val < min) min = val;
        if(val > max) max = val;
    }

    double step = (max - min) / (bins-1);

    for(const float val : bond.values_){
        int loc = int((val - min) / step);
        if(loc < 0 || loc > bins-1) cout << loc << endl;
        histogram_(loc)++;
    }
}

// use these properties to form the gaussian and check R2 value
// if it's small, go on and calculate the force constant, otherwise... use guess??
// I don't need to calculate skew or kurtosis unless they're useful for uni/multi-modal check
void BoltzmannInverter::statisticalMoments(const vector<float> &vec){
//    if(array.getDimensions() != 1) throw std::logic_error("Can't get moments of n-dimensional array");
    double sum = 0.0;
//    int n = array.getElems();
    const unsigned long n = vec.size();

    // calculate mean with first pass
//    for(int i = 0; i < n; i++) sum += array(i);
    for(int i = 0; i < n; i++) sum += vec[i];
    const double avg = sum / n;

    // calculate other stats with second pass
    double adev = 0.0, var = 0.0, skew = 0.0, kurt = 0.0, ep = 0.0;
    for(int i = 0; i < n; i++){
//        sum = array(i) - avg;
        sum = vec[i] - avg;
        ep += sum;
        adev += fabs(sum);

        double p = sum * sum;
        var += p;

        p *= sum;
        skew += p;

        p *= sum;
        kurt += p;
    }

    adev /= n;
    var = (var - ep*ep/n) / (n - 1);
    const double sdev = sqrt(var);

    if(var > 0.01){
        skew /= n * sdev*sdev*sdev;
        kurt /= (n * var*var) - 3;
    }
}