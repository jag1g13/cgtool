//
// Created by james on 10/07/15.
//

#include "light_array.h"

#include <stdio.h>
#include <stdexcept>

/** \brief Copy constructor */
template <typename T> LightArray<T>::LightArray(const LightArray &other){
    alloc(other.size_[0], other.size_[1]);
    for(int i=0; i<length_; i++){
        array_[i] = other.array_[i];
    }
}

/** \brief Assignment operator */
template <typename T> LightArray<T>& LightArray<T>::operator=(const LightArray<T> &other){
   if(size_[0]==other.size_[0] && size_[1]==other.size_[1]){
       for(int i=0; i<length_; i++){
           array_[i] = other.array_[i];
       }
   }
}

template <typename T> void LightArray<T>::alloc(const int x, const int y){
    size_[0] = x; size_[1] = y;
    length_ = x * y;
    array_ = new T[length_];
    if(array_ == nullptr) throw std::runtime_error("Could not allocate array");
}

template <typename T> T& LightArray<T>::operator()(const int x, const int y){
    return array_[x*size_[1] + y];
}

template <typename T> const T& LightArray<T>::at(const int x, const int y) const{
    return array_[x*size_[1] + y];
}

/** \brief Apply n iterations of Jacobi smoothing */
template <typename T> void LightArray<T>::smooth(const int n_iter){
    int jsw = 0;
    int isw = 0;
    for(int ipass = 0; ipass < 2 * n_iter; ipass++){
        jsw = isw;
        for(int i = 1; i < size_[0]-1; i++){
            for(int j = jsw+1; j < size_[1]-1; j += 2){
                const int loc = i*size_[0] + j;
                array_[loc] += 0.25 * (array_[loc+size_[0]] + array_[loc+1] + array_[loc-1] + array_[loc-size_[0]] - 4 * array_[loc]);
            }
            jsw = 1 - jsw;
        }
        isw = 1 - isw;
    }
}

/** \brief Print the array - format required to account for different types */
template <typename T> void LightArray<T>::print(const char *format) const{
    for(int i=0; i<size_[0]; i++){
        for(int j=0; j<size_[1]; j++){
            printf(format, array_[i*size_[0] + j]);
        }
        printf("\n");
    }
}
