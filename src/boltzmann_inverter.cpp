#include "boltzmann_inverter.h"

#include <flens/flens.cxx>

using std::cout;
using std::endl;
using std::vector;

BoltzmannInverter::BoltzmannInverter(const double temp, const int bins){
    if(bins != -1) bins_ = bins;
    temp_ = temp;
    histogram_.init(bins_);
    gaussian_.init(bins_);
    harmonic_.init(bins_);
}

void BoltzmannInverter::calculate(BondStruct &bond){
    histogram_.zero();
    gaussian_.zero();
    harmonic_.zero();
    n_ = bond.values_.size();
    statisticalMoments(bond.values_);
    bond.avg_ = mean_;
    binHistogram(bond.values_);
    bond.rsqr_ = gaussianRSquared();
    type_ = bond.type_;
    bond.forceConstant_ = invertGaussian();
}

double BoltzmannInverter::invertGaussian(){
    // R in kJ.K-1.mol-1
    const double R = 8.314 / 1000.;

    typedef flens::GeMatrix<flens::FullStorage<double>> GeMatrix;
    typedef flens::DenseVector<flens::Array<double>> DenseVector;

    GeMatrix A(bins_, 3);
    DenseVector b(bins_);

    // Setup matrices
    double x = min_ + 0.5*step_;
    double cosx;
    for(int i=1; i<=bins_; i++){
        switch(type_){
            case BondType::LENGTH:
                A(i, 1) = 1.;
                A(i, 2) = x;
                A(i, 3) = x*x;
                b(i) = -R * temp_ * log(gaussian_(i - 1) / (x * x));
                break;
            case BondType::DIHEDRAL:
            case BondType::ANGLE:
                // Angles in GROMACS are a cos^2 term
                cosx = cos(x);
                A(i, 1) = 1.;
                A(i, 2) = cosx;
                A(i, 3) = cosx*cosx;
                b(i) = -R * temp_ * log(gaussian_(i - 1) / (cosx * cosx));
                break;
//            case BondType::DIHEDRAL:
                // Dihedrals in GROMACS are k * (1 + cos(n * x - x_min))
//                cosx = 1 + cos(x);
//                A(i, 1) = 1.;
//                A(i, 2) = cosx;
//                A(i, 3) = 0;
//                b(i) = -R * T * log(gaussian_(i - 1) / (cosx));
//                break;
        };
        harmonic_(i-1) = b(i);
        x += step_;
    }

    // Solve least squares using LAPACK
    flens::lapack::ls(flens::NoTrans, A, b);
    // Sometimes gives negative force constants - fix this
    return fabs(b(3));

//    switch(type_){
//        case BondType::LENGTH:
//        case BondType::ANGLE:
//            return b(3);
//        case BondType::DIHEDRAL:
//            return b(2);
//    };
}

void BoltzmannInverter::binHistogram(const vector<double> &vec){
    step_ = (max_ - min_) / (bins_-1);
    meanBin_ = int((mean_ - min_) / step_);

    int loc = 0;
    for(const double val : vec){
        loc = int((val - min_) / step_);
        if(loc < 0 || loc > bins_-1) cout << loc << endl;
        histogram_(loc)++;
    }
}

double BoltzmannInverter::gaussianRSquared(){
    const double prefactor = 1. / (sdev_ * sqrt(2. * M_PI));
    const double postfactor = 1. / (2. * var_);

    double y_bar = 0.;
    // First pass to calculate mean and gaussian integral
    for(int i=0; i<bins_; i++){
        double x = min_ + (i + 0.5) * step_;
        double gau = prefactor * exp(-(x-mean_) * (x-mean_) * postfactor);
        gaussian_(i) = gau;
        y_bar += int(histogram_(i));
    }

    y_bar /= bins_;
    integral_ = n_ / gaussian_.sum();

    double ss_res = 0., ss_reg = 0., ss_tot = 0.;
    double sse = 0.;
    // Second pass to calculate R^2
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
    return r_sqr;
}

void BoltzmannInverter::statisticalMoments(const vector<double> &vec){
    double sum = 0.;
    // Calculate mean with first pass
    for(const double val : vec) sum += val;
    mean_ = sum / n_;
    max_ = mean_; min_ = mean_;

    double ep = 0.;
    var_ = 0.; adev_ = 0.;
    // Calculate deviations with second pass
    for(const double val : vec){
        const double dev = val - mean_;
        ep += dev;
        adev_ += fabs(dev);
        var_ += dev * dev;
        if(val < min_) min_ = val;
        if(val > max_) max_ = val;
    }

    adev_ /= n_;
    var_ = var_ / (n_ - 1);
    sdev_ = sqrt(var_);
}

void BoltzmannInverter::printGraph(Array &arr, const int scale){
    int max_num = 0;
    for(int i=0; i<bins_; i++){
        if(int(arr(i)) > max_num) max_num = int(arr(i));
    }

    // Go down rows in terminal and print marker if h_ is greater
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
