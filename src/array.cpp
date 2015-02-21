#include <stdexcept>
#include <iostream>
#include <cassert>

#include <math.h>

#include "array.h"

using std::vector;
using std::cout;
using std::endl;

ArrayFloat::ArrayFloat(){
}

//TODO turn this into a copy constructor
ArrayFloat::ArrayFloat(const int a, const int b, const int c, const bool fast){
    fast_ = fast;
    allocated_ = false;
    assert(a > 0);
    assert(b > 0);
    assert(c > 0);
    dimensions_ = 0;
    if(a > 1) dimensions_++;
    if(b > 1) dimensions_++;
    if(c > 1) dimensions_++;
    size_.reserve(3);
    size_[0] = a; size_[1] = b; size_[2] = c;
    sizex_ = a; sizey_ = b; sizez_ = c;
    elems_ = a*b*c;
    array_ = new float[elems_];
    if(array_ == NULL) throw std::runtime_error("Array alloc failed");
    allocated_ = true;
    zero();
}

void ArrayFloat::init(const int a, const int b, const int c, const bool fast){
    fast_ = fast;
    allocated_ = false;
    assert(a > 0);
    assert(b > 0);
    assert(c > 0);
    dimensions_ = 0;
    if(a > 1) dimensions_++;
    if(b > 1) dimensions_++;
    if(c > 1) dimensions_++;
    size_.reserve(3);
    size_[0] = a; size_[1] = b; size_[2] = c;
    sizex_ = a; sizey_ = b; sizez_ = c;
    elems_ = a*b*c;
    array_ = new float[elems_];
    if(array_ == NULL) throw std::runtime_error("Array alloc failed");
    allocated_ = true;
    zero();
}

void ArrayFloat::append(const vector<float> &vec){
    if(!fast_) {
        assert(allocated_);
        assert(dimensions_ == 2);
        assert(vec.size() <= size_[1]);
        assert(appendedRows_ < size_[0]);
    }
    for(int i=0; i<vec.size(); i++){
        array_[appendedRows_*sizey_ + i] = vec[i];
    }
    appendedRows_++;
}

void ArrayFloat::append(const float *vec, const int len){
    if(!fast_){
        assert(allocated_);
        assert(dimensions_ == 2);
        assert(len <= size_[1]);
        assert(appendedRows_ < size_[0]);
    }
    for(int i=0; i<len; i++){
        array_[appendedRows_*sizey_ + i] = vec[i];
    }
    appendedRows_++;
}

float& ArrayFloat::operator()(int x){
    if(!fast_) {
        assert(allocated_);
        assert(dimensions_ >= 1);
        if(x < 0) x = size_[0] + x;
        assert(x < size_[0] && x >= 0);
        /* if 2d array return ref to a row */
        if(dimensions_ == 2){
            return array_[x * size_[1]];
        }
    }
    return array_[x];
}

//float* ArrayFloat::operator()(int x){
//    if(!fast_) {
//        assert(allocated_);
//        assert(dimensions_ >= 1);
//        if(x < 0) x = size_[0] + x;
//        assert(x < size_[0] && x >= 0);
//        /* if 2d array return ref to a row */
//        if(dimensions_ == 2){
//            return array_ + (x * size_[1]);
//        }
//    }
//    return array_ + x;
//}

float& ArrayFloat::operator()(int x, int y) {
    if(fast_){
        return array_[x * sizey_ + y];
    }else{
        assert(allocated_);
        assert(dimensions_ == 2 || dimensions_ == 3);
        if(x < 0) x = size_[0] + x;
        if(y < 0) y = size_[1] + y;
        assert(x < size_[0] && x >= 0);
        assert(y < size_[1] && y >= 0);
        /* if 3d array return ref to a row */
        if(dimensions_ == 3) return array_[x * size_[1] * size_[2] + y * size_[2]];
        return array_[x * sizey_ + y];
    }
}

float& ArrayFloat::operator()(int x, int y, int z){
    if(!fast_){
        assert(allocated_);
        assert(dimensions_ == 3);
        if(x < 0) x = size_[0] + x;
        if(y < 0) y = size_[1] + y;
        if(z < 0) z = size_[2] + z;
        assert(x < size_[0] && x >= 0);
        assert(y < size_[1] && y >= 0);
        assert(z < size_[2] && z >= 0);
    }
    return array_[x * sizey_ * sizez_ + y * sizez_ + z];
}

void ArrayFloat::linspace(const int a, const int n, const float min, const float max){
    assert(allocated_);
    assert(dimensions_ == 2);
    assert(a < size_[0]);
    float *tmp = array_ + a*size_[1];
    for(int i=0; i<n; i++){
        tmp[i] = min + i*(max-min)/(n-1);
    }
}

void ArrayFloat::linspace(const int a, const int b, const int n, const float min, const float max){
    assert(allocated_);
    assert(dimensions_ == 3);
    assert(a < size_[0]);
    assert(b < size_[1]);
    float *tmp = array_ + a*size_[1]*size_[2] + b*size_[2];
    for(int i=0; i<n; i++){
        tmp[i] = min + i*(max-min)/(n-1);
    }
}

void ArrayFloat::linspace(const int n, const float min, const float max){
    assert(allocated_);
    assert(dimensions_ == 1);
    for(int i=0; i<n; i++){
        array_[i] = min + i*(max-min)/(n-1);
    }
}

void ArrayFloat::zero(){
    assert(allocated_);
    for(int i=0; i<elems_; i++){
        array_[i] = 0.f;
    }
}

void ArrayFloat::print(const int width, const int prec, const float scale){
    assert(allocated_);
    if(dimensions_ == 3) throw std::runtime_error("Not implemented");

    if(dimensions_ == 1){
        for(int i=0; i<sizex_; i++){
            printf("%*.*f", width, prec, scale*array_[i]);
        }
        cout << endl;
    }

    else if(dimensions_ == 2){
        for(int i=0; i<sizex_; i++){
            for(int j=0; j<sizey_; j++){
                printf("%*.*f", width, prec, scale*array_[i*sizey_ + j]);
            }
            cout << endl;
        }
    }
}

void ArrayFloat::free(){
    if(array_) delete[] array_;
    allocated_ = false;
    array_ = nullptr;
}

ArrayFloat::~ArrayFloat(){
    free();
}

double ArrayFloat::sum(){
    assert(allocated_);
    double sum = 0.;
    for(int i=0; i<elems_; i++) sum += array_[i];
    return sum;
}


bool operator==(const ArrayFloat &a, const ArrayFloat &b){
    assert(a.allocated_);
    assert(b.allocated_);
    assert(a.elems_ == b.elems_);
    bool equal = true;

    // if any elements are different, arrays are different
    for(int i=0; i<a.elems_; i++){
        if(a.array_[i] != b.array_[i]){
            equal = false;
            break;
        }
    }
    return equal;
}

float rms(const ArrayFloat *a, const ArrayFloat *b){
    // Both arrays need to be allocated and the same size
    assert(a->allocated_);
    assert(b->allocated_);
    assert(a->elems_ == b->elems_);
    float sum = 0.f;
    for(int i=0; i<a->elems_; i++){
        sum += (a->array_[i] - b->array_[i]) * (a->array_[i] - b->array_[i]);
    }
    return float(sqrt(sum / a->elems_));
}

float vector_rms(const vector<float> *a, const vector<float> *b){
    assert(a->size() == b->size());
    assert(a->size() != 0);
    assert(b->size() != 0);
    float sum = 0.f;
    for(int i=0; i<a->size(); i++){
        sum += ((*a)[i] - (*b)[i]) * ((*a)[i] - (*b)[i]);
    }
    return float(sqrt(sum / a->size()));
}

ArrayFloat& ArrayFloat::operator-=(const ArrayFloat &other){
    assert(elems_ == other.elems_);
    assert(dimensions_ == other.dimensions_);
    // could assert that it's the same shape, but might be useful in the future
    for(int i=0; i<elems_; i++) array_[i] -= other.array_[i];
    return (*this);
}

ArrayFloat& ArrayFloat::operator+=(const ArrayFloat &other){
    assert(elems_ == other.elems_);
    assert(dimensions_ == other.dimensions_);
    // could assert that it's the same shape, but might be useful in the future
    for(int i=0; i<elems_; i++) array_[i] += other.array_[i];
    return (*this);
}

//TODO fix this with openmp - currently shows -nan
StatsBox vector_stats(const vector<float> *a, const vector<float> *b){
    assert(a->size() == b->size());
    const int N = a->size();
    assert(a->size() != 0);
    assert(b->size() != 0);
    StatsBox result;
    float sumsqr = 0.f;
    for(int i=0; i<a->size(); i++){
        sumsqr += ((*a)[i] - (*b)[i]) * ((*a)[i] - (*b)[i]);
        result.mean_a += (*a)[i];
        result.mean_b += (*b)[i];
    }
    result.mean_a /= N;
    result.mean_b /= N;
    result.mean_diff = float(fabs(result.mean_a - result.mean_b));
    result.rms = float((sqrt(sumsqr / N)));
    result.rrms = float(fabs(result.rms / (result.mean_a)));
    return result;
}
