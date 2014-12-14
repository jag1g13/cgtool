#include <stdexcept>
#include <iostream>

#include "general.h"

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
    array_float_3d array = (float***)malloc(a * sizeof(float*));
    if(array == NULL) throw std::runtime_error("Array alloc failed");
    array[0] = (float**)malloc(a * b * sizeof(float));
    if(array[0] == NULL) throw std::runtime_error("Array alloc failed");
    array[0][0] = (float*)malloc(a * b * c * sizeof(float));
    if(array[0][0] == NULL) throw std::runtime_error("Array alloc failed");

    for(int i=0; i<b; i++){
        array[i] = array[0] + i * b * c;
        for(int j=0; j<c; j++) {
            //std::cout << "i: " << i << " j: " << j << std::endl;
            array[i][j] = array[0][0] + i * b + j * c;
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
