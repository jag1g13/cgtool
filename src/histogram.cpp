//
// Created by james on 26/05/15.
//

#include "histogram.h"

#include <cassert>
#include <stdexcept>

#include <cstdio>

Histogram::Histogram(const int size, const bool fast){
    init(size, fast);
}

Histogram::~Histogram(){
    free();
}

void Histogram::init(const int size, const bool fast){
    assert(size > 0);
    if(array_) free();

    fast_ = fast;
    allocated_ = false;

    size_ = size;
    array_ = new int[size];
    if(array_ == nullptr) throw std::runtime_error("Array alloc failed - Histogram");
    allocated_ = true;
    if(!fast_) zero();
}

void Histogram::zero(const int set){
    // Where set is zero by default - but can be specified
    for(int i=0; i<size_; i++) array_[i] = set;
}

void Histogram::free(){
    if(array_) delete[] array_;
    allocated_ = false;
    array_ = nullptr;
}

void Histogram::increment(int loc){
    if(fast_) array_[loc]++;

    assert(loc < size_);
    if(loc < 0) loc = size_ + loc;
    assert(loc >= 0);

    array_[loc]++;
}

void Histogram::decrement(int loc){
    if(fast_) array_[loc]--;

    assert(loc < size_);
    if(loc < 0) loc = size_ + loc;
    assert(loc >= 0);

    array_[loc]--;
}

int Histogram::at(int loc) const{
    if(fast_) return array_[loc];

    assert(loc < size_);
    if(loc < 0) loc = size_ + loc;
    assert(loc >= 0);

    return array_[loc];
}

void Histogram::scale(const double mult){
    for(int i=0; i<size_; i++) array_[i] *= mult;
}

void Histogram::print(const int width) const{
    assert(allocated_);

    for(int i=0; i<size_; i++) printf("%*d", width, array_[i]);
    printf("\n");
}

void Histogram::printGraph(const int scale) const{
    int max_num = 0;
    for(int i=0; i<size_; i++){
        if(array_[i] > max_num) max_num = array_[i];
    }

    const double bar_scale = static_cast<double>(scale) / max_num;

    // Go down rows in terminal and print marker if array_ is greater
    for(int i=scale; i>0; i--){
        printf("%5.3f|", min_);
        for(int j=0; j<size_; j++){
            if(array_[j]*bar_scale >= i){
                printf("#");
            }else{
                printf(" ");
            }
        }
        printf("|%5.3f\n", max_);
    }
}
