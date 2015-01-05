#include <stdexcept>
#include <iostream>
#include <cassert>
//#include <exception>

#include "arrays.h"

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
    cout << size_[0] << "x" << size_[1] << "x" << size_[2] << endl;
    elems_ = a*b*c;
//    for(int i : size_) elems_ *= i;
    array_ = new float[elems_];
    if(array_ == NULL) throw std::runtime_error("Array alloc failed");
    allocated_ = true;
    cout << dimensions_ << "d array\t";
//    for(int i : size_) cout << i << "x";
    cout << elems_ << " elements" << endl;
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
    cout << size_[0] << "x" << size_[1] << "x" << size_[2] << endl;
    elems_ = a*b*c;
//    for(int i : size_) elems_ *= i;
    array_ = new float[elems_];
    if(array_ == NULL) throw std::runtime_error("Array alloc failed");
    allocated_ = true;
    cout << dimensions_ << "d array\t";
//    for(int i : size_) cout << i << "x";
    cout << elems_ << " elements" << endl;
}

void ArrayFloat::append(vector<float> vec){
    if(!fast_) {
        assert(dimensions_ == 2);
        assert(vec.size() <= size_[1]);
        assert(appendedRows_ < size_[0]);
    }
    for(int i=0; i<vec.size(); i++){
        array_[appendedRows_*sizey_ + i] = vec[i];
    }
    appendedRows_++;
}

float& ArrayFloat::operator()(int x){
    if(!fast_) {
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

float& ArrayFloat::operator()(int x, int y) {
    //return array_[x * size_[1] + y];
    if(fast_){
        return array_[x * sizey_ + y];
    }else{
        assert(dimensions_ == 2 || dimensions_ == 3);
//        if (x < 0) x = size_[0] + x;
//        if (y < 0) y = size_[1] + y;
        assert(x < size_[0] && x >= 0);
        assert(y < size_[1] && y >= 0);
        /* if 3d array return ref to a row */
        if (dimensions_ == 3) {
            return array_[x * size_[1] * size_[2] + y * size_[2]];
        }
        return array_[x * sizey_ + y];
    }
}

float& ArrayFloat::operator()(int x, int y, int z){
    if(!fast_) {
        assert(dimensions_ == 3);
        if(x < 0) x = size_[0] + x;
        if(y < 0) y = size_[1] + y;
        if(z < 0) z = size_[2] + z;
        assert(x < size_[0] && x >= 0);
        assert(y < size_[1] && y >= 0);
        assert(z < size_[2] && z >= 0);
    }
    //return array_[x * size_[1] * size_[2] + y * size_[2] + z];
    return array_[x * sizey_ * sizez_ + y * sizez_ + z];
}

void ArrayFloat::linspace(const int a, const int n, const float min, const float max){
    assert(dimensions_ == 2);
    assert(a < size_[0]);
    float *tmp = array_ + a*size_[1];
    for(int i=0; i<n; i++){
        tmp[i] = min + i*(max-min)/(n-1);
    }
}

void ArrayFloat::linspace(const int a, const int b, const int n, const float min, const float max){
    assert(dimensions_ == 3);
    assert(a < size_[0]);
    assert(b < size_[1]);
    float *tmp = array_ + a*size_[1]*size_[2] + b*size_[2];
    for(int i=0; i<n; i++){
        tmp[i] = min + i*(max-min)/(n-1);
    }
}

void ArrayFloat::linspace(const int n, const float min, const float max){
    assert(dimensions_ == 1);
    for(int i=0; i<n; i++){
        array_[i] = min + i*(max-min)/(n-1);
    }
}

void ArrayFloat::zero(){
    for(int i=0; i<elems_; i++){
        array_[i] = 0.;
    }
}

//ArrayFloat::~ArrayFloat(){
//    if(array_!=nullptr){
//        cout << "Freeing " << this << endl;
//        //free(array_);
//        delete[] array_;
//        array_ = nullptr;
//    }
//}

void linspace_1d(array_float_1d array, const float min, const float max, const int steps){
    for(int i=0; i<steps; i++){
        array[i] = min + steps*(max-min);
    }
}

