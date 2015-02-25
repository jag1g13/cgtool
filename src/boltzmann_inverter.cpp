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
    const double R = 8.314;
    const int T = 310;

    ArrayFloat harmonic(bins_);

    double x = min_ + 0.5*step_;
    for(int i=0; i<bins_; i++){
        harmonic(i) = -R * T * log(gaussian_(i) / (x*x));
        x += step_;
    }

    // equation y = kx(x - mean)

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

    for(const double val : bond.values_){
        int loc = int((val - min_) / step_);
        if(loc < 0 || loc > bins-1) cout << loc << endl;
        histogram_(loc)++;
    }
}

double BoltzmannInverter::gaussianRSquared(){
    const double prefactor = 1. / (sdev_ * sqrt(M_2_PI));
    const double postfactor = 1. / (2 * var_);
    double y_bar = 0.;

    // first pass to calculate mean and gaussian integral
    for(int i=0; i<bins_; i++){
        double x = min_ + (i + 0.5) * step_;
        double gau = prefactor * exp(-(x-mean_) * (x-mean_) * postfactor);
        gaussian_(i) = gau;
        y_bar += int(histogram_(i));
    }

    y_bar /= bins_;
    const double gau_scale = n_ / gaussian_.sum();
    double ss_res = 0., ss_reg = 0., ss_tot = 0.;
    double sse = 0.;

    // second pass to calculate R^2
    for(int i=0; i<bins_; i++){
        int actual = int(histogram_(i));
        double gau = gaussian_(i) * gau_scale;
//        ss_reg += (gau - y_bar) * (gau - y_bar);
        ss_res += (gau - actual) * (gau - actual);
        ss_tot += (actual - y_bar) * (actual - y_bar);
//        sse += (actual - gau) * (actual - gau);
    }
    const double r_sqr = 1 - ss_res / ss_tot;
//    const double r_sqr2 = ss_reg / ss_tot;
//    sse = log(sse / n_);
//    printf("%8.3f%8.3f%12.3f\n", r_sqr, r_sqr2, sse);
    printf("%8.3f\n", r_sqr);
    return r_sqr;
}

void BoltzmannInverter::statisticalMoments(const vector<double> &vec){
    double sum = 0.0;
    const unsigned long n = vec.size();

    // calculate mean with first pass
    for(int i = 0; i < n; i++) sum += vec[i];
    mean_ = sum / n;

    // calculate deviations with second pass
    double ep = 0.0;
    var_ = 0.0; adev_ = 0.0;
    for(int i = 0; i < n; i++){
        double dev = vec[i] - mean_;
        ep += dev;
        adev_ += fabs(dev);
        var_ += dev * dev;
    }

    adev_ /= n;
    var_ = (var_ - ep*ep/n) / (n - 1);
    sdev_ = sqrt(var_);
}
