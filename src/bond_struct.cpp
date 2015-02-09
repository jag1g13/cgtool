#include "bond_struct.h"

#include <iostream>

using std::cout;
using std::endl;

int BondStruct::totalNum_ = 0;

BondStruct::BondStruct(const int size){
    // keep track of serial number
    num_ = totalNum_;
    totalNum_++;
    atomNames_.resize(size);
    atomNums_.resize(size);
}

BondStruct::BondStruct(const BondStruct &other){
    unsigned long size = other.atomNames_.size();
    atomNames_.resize(size);
    atomNums_.resize(size);
    for(int i = 0; i < size; i++){
        atomNames_[i] = other.atomNames_[i];
        atomNums_[i] = other.atomNums_[i];
    }
}

void BondStruct::calcAvg(){
    for(const float &val : values_){
        avg_ += val;
    }
    avg_ /= values_.size();
}
