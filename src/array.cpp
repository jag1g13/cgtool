#include <stdexcept>
#include <iostream>
#include <cassert>

#include <math.h>

#include "array.h"

using std::vector;
using std::cout;
using std::endl;
using std::string;

//TODO turn this into a copy constructor
Array::Array(const int a, const int b, const int c, const bool fast){
    init(a, b, c, fast);
}

void Array::init(const int a, const int b, const int c, const bool fast){
    if(array_) free();

    fast_ = fast;
    allocated_ = false;

    assert(a > 0);
    assert(b > 0);
    assert(c > 0);

    dimensions_ = 0;
    if(a > 1) dimensions_++;
    if(b > 1) dimensions_++;
    if(c > 1) dimensions_++;

    size_[0] = a; size_[1] = b; size_[2] = c;
    elems_ = a*b*c;
    array_ = new double[elems_];

    if(array_ == NULL) throw std::runtime_error("Array alloc failed");
    allocated_ = true;
    if(!fast) zero();
}

void Array::append(const vector<double> &vec){
    if(!fast_) {
        assert(allocated_);
        assert(dimensions_ == 2);
        assert(vec.size() <= size_[1]);
        assert(appendedRows_ < size_[0]);
    }
    for(int i=0; i<vec.size(); i++){
        array_[appendedRows_*size_[1] + i] = vec[i];
    }
    appendedRows_++;
}

void Array::append(const double *vec, int len){
    if(len == -1) len = size_[1];
    if(!fast_) {
        assert(allocated_);
        assert(dimensions_ == 2);
        assert(len <= size_[1]);
        assert(appendedRows_ < size_[0]);
    }
    for(int i=0; i<len; i++){
        array_[appendedRows_*size_[1] + i] = vec[i];
    }
    appendedRows_++;
}

double& Array::operator()(int x){
    if(fast_) return array_[x];

    assert(allocated_);
    assert(dimensions_ >= 1);
    if(x < 0) x = size_[0] + x;
    assert(x < size_[0] && x >= 0);
    /* if 2d array return ref to a row */
    if(dimensions_ == 2) return array_[x * size_[1]];
    /* if 3d array return ref to a plane */
    if(dimensions_ == 3) return array_[x * size_[1] * size_[2]];
    return array_[x];
}

double& Array::operator()(int x, int y){
    if(fast_) return array_[x * size_[1] + y];

    assert(allocated_);
    assert(dimensions_ >= 2);
    if(x < 0) x = size_[0] + x;
    if(y < 0) y = size_[1] + y;
    assert(x < size_[0] && x >= 0);
    assert(y < size_[1] && y >= 0);
    /* if 3d array return ref to a row */
    if(dimensions_ == 3) return array_[x * size_[1] * size_[2] + y * size_[2]];
    return array_[x * size_[1] + y];
}

double& Array::at(int x, int y) {
    assert(allocated_);
    assert(dimensions_ == 2);
    x = x % size_[0];
    y = y % size_[1];
    while(x < 0) x += size_[0];
    while(y < 0) y += size_[1];
    return array_[x * size_[1] + y];
}


double& Array::operator()(int x, int y, int z){
    if(fast_) return array_[x * size_[1] * size_[2] + y * size_[2] + z];

    assert(allocated_);
    assert(dimensions_ == 3);
    if(x < 0) x = size_[0] + x;
    if(y < 0) y = size_[1] + y;
    if(z < 0) z = size_[2] + z;
    assert(x < size_[0] && x >= 0);
    assert(y < size_[1] && y >= 0);
    assert(z < size_[2] && z >= 0);
    return array_[x * size_[1] * size_[2] + y * size_[2] + z];
}

void Array::linspace(const int n, const double min, const double max){
    assert(allocated_);
    assert(dimensions_ == 1);
    assert(n <= size_[0]);
    const double step = (max - min) / (n - 1);

    for(int i=0; i<n; i++) array_[i] = min + i*step;
}


void Array::linspace(const int a, const int n, const double min, const double max){
    assert(allocated_);
    assert(dimensions_ == 2);
    assert(a < size_[0]);
    assert(n <= size_[1]);
    double *tmp = array_ + a*size_[1];
    const double step = (max - min) / (n - 1);

    for(int i=0; i<n; i++) tmp[i] = min + i*step;
}

void Array::linspace(const int a, const int b, const int n, const double min, const double max){
    assert(allocated_);
    assert(dimensions_ == 3);
    assert(a < size_[0]);
    assert(b < size_[1]);
    double *tmp = array_ + a*size_[1]*size_[2] + b*size_[2];
    const double step = (max - min) / (n - 1);

    for(int i=0; i<n; i++) tmp[i] = min + i*step;
}

void Array::zero(){
    assert(allocated_);
    for(int i=0; i<elems_; i++) array_[i] = 0.;
}

void Array::print(const int width, const int prec, const double scale){
    assert(allocated_);
    // Print if 1d or 2d, otherwise ignore

    if(dimensions_ == 1){
        for(int i=0; i<size_[0]; i++){
            printf("%*.*f", width, prec, scale * array_[i]);
        }
        cout << endl;
    }

    else if(dimensions_ == 2){
        for(int i=0; i<size_[0]; i++){
            for(int j=0; j<size_[1]; j++){
                printf("%*.*f", width, prec, scale*array_[i*size_[1] + j]);
            }
            cout << endl;
        }
    }
}

void Array::printCSV(const std::string &filename, const int remove_border){
    const int r = remove_border;
    assert(r >= 0);
    const string file = filename + ".dat";

    // Backup using small_functions.h
//    backup_old_file(file);

    FILE *f = fopen(file.c_str(), "a");

    if(dimensions_ == 1){
        for(int i=r; i < size_[0]-r; i++) fprintf(f, "%8.3f\n", array_[i]);

    }else if(dimensions_ == 2){
        for(int i=r; i < size_[0]-r; i++){
            if(dimensions_ == 2){}
            for(int j=r; j < size_[1]-r; j++){
                fprintf(f, "%8.3f", array_[i*size_[0] + j]);
            }
            fprintf(f, "\n");
        }
    }
    // If 3d do nothing

    fclose(f);
}

void Array::free(){
    // Check for nullptr - has it been (de)allocated yet
    if(array_) delete[] array_;
    allocated_ = false;
    array_ = nullptr;
}

Array::~Array(){
    free();
}

double Array::sum(){
    assert(allocated_);
    double sum = 0.;
    for(int i=0; i<elems_; i++) sum += array_[i];
    return sum;
}

double Array::mean(){
    assert(allocated_);
    return sum() / elems_;
}

double rmsd(const Array &a, const Array &b){
    // Both arrays need to be allocated and the same size - not shape
    assert(a.allocated_);
    assert(b.allocated_);
    assert(a.elems_ == b.elems_);
    double sum = 0.;
    for(int i=0; i<a.elems_; i++){
        sum += (a.array_[i] - b.array_[i]) * (a.array_[i] - b.array_[i]);
    }
    return sqrt(sum / a.elems_);
}

bool operator==(const Array &a, const Array &b){
    // Both arrays need to be allocated
    assert(a.allocated_);
    assert(b.allocated_);
    if(a.elems_ != b.elems_) return false;
    bool equal = true;

    // If any elements are different, arrays are different
    for(int i=0; i<a.elems_; i++){
        if(a.array_[i] != b.array_[i]){
            equal = false;
            break;
        }
    }
    return equal;
}

Array& Array::operator-=(const Array &other){
    assert(elems_ == other.elems_);
    assert(dimensions_ == other.dimensions_);
    for(int i=0; i<elems_; i++) array_[i] -= other.array_[i];
    return (*this);
}

Array& Array::operator+=(const Array &other){
    assert(elems_ == other.elems_);
    assert(dimensions_ == other.dimensions_);
    for(int i=0; i<elems_; i++) array_[i] += other.array_[i];
    return (*this);
}

Array& Array::operator*=(const double mult){
    for(int i=0; i<elems_; i++) array_[i] *= mult;
    return (*this);
}

Array& Array::operator/=(const double div){
    const double mult = 1. / div;
    for(int i=0; i<elems_; i++) array_[i] *= mult;
    return (*this);
}

Array& Array::operator+=(const double add){
    for(int i=0; i<elems_; i++) array_[i] += add;
    return (*this);
}

Array& Array::operator-=(const double sub){
    for(int i=0; i<elems_; i++) array_[i] -= sub;
    return (*this);
}

void Array::elementMultiply(const Array &other){
    assert(elems_ == other.elems_);
    assert(dimensions_ == other.dimensions_);
    for(int i=0; i<elems_; i++) array_[i] *= other.array_[i];
}

void Array::elementDivide(const Array &other){
    assert(elems_ == other.elems_);
    assert(dimensions_ == other.dimensions_);

    for(int i=0; i<elems_; i++) array_[i] /= other.array_[i];
}

void Array::replaceNaN(){
    vector<int> nans;
    for(int i=0; i<elems_; i++){
        // If NaN - uses the fact that NaN != NaN
        if(array_[i] != array_[i]){
            array_[i] = 0.;
            nans.push_back(i);
        }
    }

    const double m = mean();
    for(int &i : nans) array_[i] = m;
}

void Array::smooth(const int n_iter){
    int jsw = 0;
    int isw = 0;
    for(int ipass=0; ipass < 2 * n_iter; ipass++){
        jsw = isw;
        for(int i=0; i < size_[0]; i++){
            for(int j=jsw; j < size_[1]; j+=2){
                at(i, j) += 0.25 * (at(i-1, j) + at(i,j+1) + at(i,j-1) + at(i+1, j) - 4*at(i, j));
            }
            jsw = 1 - jsw;
        }
        isw = 1 - isw;
    }
}

void Array::interpolateZeros(){
    for(int i=0; i<size_[0]; i++){
        for(int j=0; j<size_[1]; j++){
            if(at(i, j) == 0.) at(i, j) = (at(i,j-1) + at(i-1, j) +
                                           at(i,j+1) + at(i+1, j)) / 4.;
        }
    }
}