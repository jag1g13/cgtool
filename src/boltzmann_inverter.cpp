#include "boltzmann_inverter.h"

#include <iostream>

#include <math.h>

#include "array.h"

using std::cout;
using std::endl;

void BoltzmannInverter::invertGaussian(){
    /* line from Python version
    y_inv = -R * T * np.log(y_fit / (x_fit*x_fit))
     */
}

void BoltzmannInverter::binHistogram(const BondStruct &bond, const int bins){
    // make sure you know where the point estimates for the
    // histogram for the residual calculation are coming from
    //
    // look up constant variance test to see if you need more data
    max_ = bond.avg_; min_ = bond.avg_;
    bins_ = bins;
    histogram_.init(bins);
    gaussian_.init(bins);
    n_ = bond.values_.size();

    for(const double val : bond.values_){
        if(val < min_) min_ = val;
        if(val > max_) max_ = val;
    }

    step_ = (max_ - min_) / (bins-1);
//    printf("%8.3f%8.3f%8.5f\n", min_, max_, step_);

    for(const float val : bond.values_){
        int loc = int((val - min_) / step_);
        if(loc < 0 || loc > bins-1) cout << loc << endl;
        histogram_(loc)++;
    }
}

double BoltzmannInverter::gaussianRSquared(){
    const double root2pi = sqrt(M_2_PI);
    double y_bar = 0.;

    // first pass to calculate mean and gaussian integral
    for(int i=0; i<bins_; i++){
        double x = min_ + i * step_;
        double gau = (1 / (sdev_ * root2pi)) * exp(-(x - avg_) * (x - avg_) / (2 * sdev_ * sdev_));
        gaussian_(i) = gau;
        y_bar += int(histogram_(i));
    }

    y_bar /= bins_;
    double gau_scale = n_ / gaussian_.sum();
    double ss_res = 0., ss_reg = 0., ss_tot = 0.;

    // second pass to calculate R^2
    for(int i=0; i<bins_; i++){
        int actual = int(histogram_(i));
        double gau = gaussian_(i) * gau_scale;
        ss_res += (actual - gau) * (actual - gau);
        ss_tot += (actual - y_bar) * (actual - y_bar);
    }
    //const double r_sqr = 1 - ss_res / ss_tot;
    const double r_sqr = ss_res / ss_tot;
    printf("%8.3f\n", r_sqr);
    return r_sqr;
}

// use these properties to form the gaussian and check R2 value
// if it's small, go on and calculate the force constant, otherwise... use guess??
// I don't need to calculate skew or kurtosis unless they're useful for uni/multi-modal check
void BoltzmannInverter::statisticalMoments(const vector<double> &vec){
//    if(array.getDimensions() != 1) throw std::logic_error("Can't get moments of n-dimensional array");
    double sum = 0.0;
//    int n = array.getElems();
    const unsigned long n = vec.size();

    // calculate mean with first pass
//    for(int i = 0; i < n; i++) sum += array(i);
    for(int i = 0; i < n; i++) sum += vec[i];
    avg_ = sum / n;

    // calculate other stats with second pass
    double var = 0.0, ep = 0.0;
    adev_ = 0.0; skew_ = 0.0; kurt_ = 0.0;
    for(int i = 0; i < n; i++){
        sum = vec[i] - avg_;
        ep += sum;
        adev_ += fabs(sum);

        double p = sum * sum;
        var += p;

        p *= sum;
        skew_ += p;

        p *= sum;
        kurt_ += p;
    }

    adev_ /= n;
    var = (var - ep*ep/n) / (n - 1);
    sdev_ = sqrt(var);

    if(var > 0.01){
        skew_ /= n * sdev_*var;
        kurt_ /= (n * var*var) - 3;
    }

//    printf("%12.5f%12.5f%12.5f%12.5f\n", avg_, sdev_, skew_, kurt_);
}
