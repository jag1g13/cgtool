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

void BondStruct::calcAvg(){
    for(const float &val : values_){
        avg_ += val;
    }
    avg_ /= values_.size();
}
