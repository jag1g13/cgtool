#include <stdexcept>
#include <iostream>
#include <cassert>

#include "arrays.h"

using std::vector;
using std::cout;
using std::endl;

Array::Array(vector<int> size){
    dimensions_ = size.size();
    size_ = size;
}

ArrayFloat::ArrayFloat(){
}

ArrayFloat::ArrayFloat(const vector<int> size, const bool fast){
    fast_ = fast;
    allocated_ = false;
    dimensions_ = size.size();
    size_ = size;
    elems_ = 1;
    assert(dimensions_ <= 3);
    for(int i : size_) elems_ *= i;
    //array_ = (float*)malloc(elems_ * sizeof(float));
    array_ = new float[elems_];
    if(array_ == NULL) throw std::runtime_error("Array alloc failed");
    allocated_ = true;
    cout << dimensions_ << "d array\t";
    for(int i : size_) cout << i << " ";
    //cout << size_[0] << "x" << size_[1] << "x" << size_[2] << "\t";
    cout << elems_ << " elements" << endl;
}

float& ArrayFloat::operator()(const int x){
    if(!fast_) {
        assert(dimensions_ == 1);
        assert(x < size_[0] && x >= 0);
    }
    return array_[x];
}

float& ArrayFloat::operator()(const int x, const int y){
    if(!fast_) {
        assert(dimensions_ == 2);
        assert(x < size_[0] && x >= 0);
        assert(y < size_[1] && y >= 0);
    }
    return array_[x * size_[1]+ y];
}

float& ArrayFloat::operator()(const int x, const int y, const int z){
    if(!fast_) {
        assert(dimensions_ == 3);
        assert(x < size_[0] && x >= 0);
        assert(y < size_[1] && y >= 0);
        assert(z < size_[2] && z >= 0);
    }
    return array_[x * size_[1] * size_[2] + y * size_[2] + z];
}

void ArrayFloat::zero(){
    for(int i=0; i<elems_; i++){
        array_[i] = 0.;
    }
}

ArrayFloat::~ArrayFloat(){
    if(allocated_){
        cout << "Freeing " << this << endl;
        //free(array_);
        //delete[] array_;
    }
    allocated_ = false;
    //array_ = nullptr;
}

//unsigned long operator ArrayFloat::()(int i, int j=0, int k=0) const{
//}
//unsigned long & operator ArrayFloat::()(int i, int j=0, int k=0){
//}

array_float_1d alloc_float_3d_flat(const int a, const int b, const int c){
    array_float_1d array = (float*)malloc(a * b * c * sizeof(float));
    if(array == NULL) throw std::runtime_error("Array alloc failed");
    return array;
}

array_float_1d alloc_float_2d_flat(const int a, const int b){
    array_float_1d array = (float*)malloc(a * b * sizeof(float));
    if(array == NULL) throw std::runtime_error("Array alloc failed");
    return array;
}

array_float_2d alloc_float_2d(const int a, const int b){
    array_float_2d array = (float**)malloc(a * sizeof(float*));
    if(array == NULL) throw std::runtime_error("Array alloc failed");
    array[0] = (float*)malloc(a * b * sizeof(float));
    if(array[0] == NULL) throw std::runtime_error("Array alloc failed");

    for(int i=0; i<a; i++){
        array[i] = array[0] + i*b;
    }
    return array;
}

array_float_3d alloc_float_3d(const int a, const int b, const int c){
    //std::cout << "a:" << a << " b:" << b << " c:" << c << std::endl;
    array_float_3d array = (float***)malloc(a * sizeof(float**));
    if(array == NULL) throw std::runtime_error("Array alloc failed");
    array[0] = (float**)malloc(a * b * sizeof(float*));
    if(array[0] == NULL) throw std::runtime_error("Array alloc failed");
    array[0][0] = (float*)malloc(a * b * c * sizeof(float));
    if(array[0][0] == NULL) throw std::runtime_error("Array alloc failed");

    for(int i=0; i<a; i++){
        //array[i] = array[0] + i*b;
        array[i] = array[0] + i;
        for(int j=0; j<b; j++) {
            //std::cout << "i: " << i << " j: " << j << std::endl;
            array[i][j] = array[i][0] + j*c;
        }
    }
    return array;
}

void linspace_1d(array_float_1d array, const float min, const float max, const int steps){
    for(int i=0; i<steps; i++){
        array[i] = min + steps*(max-min);
    }
}

void zero_float_2d(array_float_2d array, const int a, const int b){
    float *tmp = array[0];
    for(int i=0; i<a*b; i++){
        tmp[i] = 0.f;
    }
}

void zero_float_3d(array_float_3d array, const int a, const int b, const int c){
    float *tmp = array[0][0];
    for(int i=0; i<a*b*c; i++){
        tmp[i] = 0.f;
    }
}
