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

    gaussian_.print(8, 2);
    double x = min_ + 0.5*step_;
    for(int i=0; i<bins_; i++){
        harmonic(i) = -R * T * log(gaussian_(i) / (x*x));
        x += step_;
    }
    harmonic.print(12, 2);
    cout << endl;

    // equation y = kx(x - mean)
    // may as well fit all coefficients
    // replace with lapack
    // A = (xT . x)^-1 . xT . y
    // where y is harmonic()
    // replace code here with matrix multiply in ArrayFloat

    // input arrays
    ArrayFloat X(bins_, 3);
    for(int i=0; i<bins_; i++){
        X(i, 0) = 1;
        X(i, 1) = (min_ + (i + 0.5) * step_);
        X(i, 2) = X(i, 1) * X(i, 1);
    }
//    X(0, 0) = 1; X(0,1) = -1; X(0, 2) = 1;
//    X(1, 0) = 1; X(1,1) = 0; X(1, 2) = 0;
//    X(2, 0) = 1; X(2,1) = 1; X(2, 2) = 1;
//    X(3, 0) = 1; X(3,1) = 2; X(3, 2) = 4;
    X.print(10, 3);
    cout << endl;

//    ArrayFloat Y(bins_);
//    Y(0) = 1; Y(1) = 0; Y(2) = 1; Y(3) = 4;

    // output array
    ArrayFloat A(3);

    // tmp arrays
    ArrayFloat XT(3, bins_);
    ArrayFloat XTX(3, 3);
    ArrayFloat XTXInv(3, 3);

    for(int i=0; i<3; i++){
        for(int j=0; j<bins_; j++){
            XT(i, j) = X(j, i);
        }
    }
//    XT.print();
//    cout << endl;

    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            for(int k=0; k<bins_; k++){
                XTX(i, j) += XT(j, k) * X(k, i);
            }
        }
    }
//    XTX.print();
//    cout << endl;

    double det = 0.;
    det += XTX(0, 0) * (XTX(1, 1)*XTX(2, 2) - XTX(2, 1)*XTX(1, 2));
    det -= XTX(1, 0) * (XTX(0, 1)*XTX(2, 2) - XTX(2, 1)*XTX(0, 2));
    det += XTX(2, 0) * (XTX(0, 1)*XTX(1, 2) - XTX(1, 1)*XTX(0, 2));
    cout << det << endl;
    cout << endl;
    double invdet = 1. / det;

    XTXInv(0,0) =  (XTX(1,1)*XTX(2,2)-XTX(2,1)*XTX(1,2))*invdet;
    XTXInv(0,1) = -(XTX(0,1)*XTX(2,2)-XTX(0,2)*XTX(2,1))*invdet;
    XTXInv(0,2) =  (XTX(0,1)*XTX(1,2)-XTX(0,2)*XTX(1,1))*invdet;
    XTXInv(1,0) = -(XTX(1,0)*XTX(2,2)-XTX(1,2)*XTX(2,0))*invdet;
    XTXInv(1,1) =  (XTX(0,0)*XTX(2,2)-XTX(0,2)*XTX(2,0))*invdet;
    XTXInv(1,2) = -(XTX(0,0)*XTX(1,2)-XTX(1,0)*XTX(0,2))*invdet;
    XTXInv(2,0) =  (XTX(1,0)*XTX(2,1)-XTX(2,0)*XTX(1,1))*invdet;
    XTXInv(2,1) = -(XTX(0,0)*XTX(2,1)-XTX(2,0)*XTX(0,1))*invdet;
    XTXInv(2,2) =  (XTX(0,0)*XTX(1,1)-XTX(1,0)*XTX(0,1))*invdet;
//    XTXInv.print();
//    cout << endl;

    ArrayFloat XTXInvXT(3, bins_);
    for(int i=0; i<3; i++){
        for(int j=0; j<bins_; j++){
            for(int k=0; k<3; k++){
                XTXInvXT(i, j) += XTXInv(i, k) * XT(k, j);
            }
        }
    }
//    XTXInvXT.print();
//    cout << endl;

    for(int i=0; i<3; i++){
        for(int k=0; k<bins_; k++){
            A(i) += XTXInvXT(i, k) * harmonic(k);
        }
    }
    A.print(12, 3);
    cout << endl;
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
        gaussian_(i) *= gau_scale;
        const double gau = gaussian_(i);
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
