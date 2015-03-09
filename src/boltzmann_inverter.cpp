#include "boltzmann_inverter.h"

#include <flens/flens.cxx>

using std::cout;
using std::endl;
using std::vector;

BoltzmannInverter::BoltzmannInverter(BondStruct &bond, const bool isAngle){
    statisticalMoments(bond.values_);
    binHistogram(bond, 55);
    bond.rsqr_ = gaussianRSquared();
//    bool isAngle = (bond.type_ != BondType::LENGTH);
    bond.forceConstant_ = invertGaussian(isAngle);
    printGraph(histogram_);
    printGraph(gaussian_);
//    printGraph(harmonic_);
//    printGraph(harmonicFit_);
}

double BoltzmannInverter::invertGaussian(const bool isAngle){
    /* line from Python version
    y_inv = -R * T * np.log(y_fit / (x_fit*x_fit))
     */
    // R in kJ.K-1.mol-1
    const double R = 8.314 / 1000.;
    const int T = 310;

    typedef flens::GeMatrix<flens::FullStorage<double> >   GeMatrix;
    typedef flens::DenseVector<flens::Array<double> >      DenseVector;
    typedef typename flens::DenseVector<flens::Array<double>>::IndexType  IndexType;

    GeMatrix     A(bins_, 3);
    DenseVector  b(bins_);

    // Setup matrices
    harmonic_.init(bins_);
    double x = min_ + 0.5*step_;
    for(int i=1; i<=bins_; i++){
        const double cosx = cos(x);
        if(isAngle){
            b(i) = -R * T * log(gaussian_(i - 1) / (cosx * cosx));
        }else{
            b(i) = -R * T * log(gaussian_(i - 1) / (x * x));
        }
        harmonic_(i-1) = b(i);
        // Angles in MARTINI are a cos^2 term
        if(isAngle){
            A(i, 1) = 1;
            A(i, 2) = cosx;
            A(i, 3) = cosx*cosx;
        }else{
            A(i, 1) = 1;
            A(i, 2) = x;
            A(i, 3) = x*x;
        }
        x += step_;
    }

    // Solve least squares using LAPACK
    flens::lapack::ls(flens::NoTrans, A, b);
    const flens::Underscore<IndexType>  _;
    auto X = b(_(1,3));
//    cout << "X = " << X << endl;
    return X(3);
}

void BoltzmannInverter::binHistogram(const BondStruct &bond, const int bins){
    // make sure you know where the point estimates for the
    // histogram for the residual calculation are coming from
    //
    // look up constant variance test to see if you need more data
    bins_ = bins;
    histogram_.init(bins);
    gaussian_.init(bins);
    n_ = bond.values_.size();


    step_ = (max_ - min_) / (bins-1);
//    printf("%8.3f%8.3f%8.5f\n", min_, max_, step_);
    meanBin_ = int((mean_ - min_) / step_);

    for(const double val : bond.values_){
        int loc = int((val - min_) / step_);
        if(loc < 0 || loc > bins-1) cout << loc << endl;
        histogram_(loc)++;
    }
}

double BoltzmannInverter::gaussianRSquared(){
//    cout << "gaussianRSquared" << endl;
    const double prefactor = 1. / (sdev_ * sqrt(2. * M_PI));
    const double postfactor = 1. / (2. * var_);
//    printf("%12.3f%12.3f%12.3f\n", amplitude_, prefactor, postfactor);
    double y_bar = 0.;

    // first pass to calculate mean and gaussian integral
    for(int i=0; i<bins_; i++){
        double x = min_ + (i + 0.5) * step_;
        double gau = prefactor * exp(-(x-mean_) * (x-mean_) * postfactor);
        gaussian_(i) = gau;
        y_bar += int(histogram_(i));
    }
//    gaussian_.print();

    y_bar /= bins_;
    integral_ = n_ / gaussian_.sum();
//    cout << integral_ << endl;
    double ss_res = 0., ss_reg = 0., ss_tot = 0.;
    double sse = 0.;

    // second pass to calculate R^2
    for(int i=0; i<bins_; i++){
        int actual = int(histogram_(i));
        gaussian_(i) *= integral_;
        const double gau = gaussian_(i);
        ss_reg += (gau - y_bar) * (gau - y_bar);
        ss_res += (gau - actual) * (gau - actual);
        ss_tot += (actual - y_bar) * (actual - y_bar);
        sse += (actual - gau) * (actual - gau);
    }
    const double r_sqr = 1 - ss_res / ss_tot;
    const double r_sqr2 = ss_reg / ss_tot;
    sse = log(sse / n_);
//    printf("%8.3f%8.3f%12.3f\n", r_sqr, r_sqr2, sse);
//    printf("%8.3f\n", r_sqr);
    return r_sqr;
}

void BoltzmannInverter::statisticalMoments(const vector<double> &vec){
    double sum = 0.0;
    const unsigned long n = vec.size();

    // calculate mean with first pass
    for(int i = 0; i < n; i++) sum += vec[i];
    mean_ = sum / n;
    max_ = mean_; min_ = mean_;

    // calculate deviations with second pass
    double ep = 0.0;
    var_ = 0.0; adev_ = 0.0;
    for(int i = 0; i < n; i++){
        double dev = vec[i] - mean_;
        ep += dev;
        adev_ += fabs(dev);
        var_ += dev * dev;
        if(vec[i] < min_) min_ = vec[i];
        if(vec[i] > max_) max_ = vec[i];
    }

    adev_ /= n;
    //TODO what's going on here - why does n/2 give a better fit
    var_ = var_ / (n - 1);
    sdev_ = sqrt(var_);
//    printf("%12.5f%12.5f%12.5f%12.5f\n", mean_, var_, sdev_, max_-min_);
}

void BoltzmannInverter::printGraph(Array &arr, const int scale){
    int max_num = 0;
    for(int i=0; i<bins_; i++){
        if(int(arr(i)) > max_num) max_num = int(arr(i));
    }
    // go down rows in terminal and print marker if h_ is greater
    for(int i=scale; i>0; i--){
        printf("%5.3f|", min_);
        for(int j=0; j<bins_; j++){
            if(arr(j)*scale/max_num >= i){
                printf("#");
            }else{
                printf(" ");
            }
        }
        printf("|%5.3f\n", max_);
    }
};
