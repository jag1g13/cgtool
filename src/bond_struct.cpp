#include "bond_struct.h"

#include <iostream>

using std::cout;
using std::endl;

BondStruct::BondStruct(const int size){
    atomNames_.resize(size);
    atomNums_.resize(size);
}

void BondStruct::calcAvg(){
    for(const float &val : values_){
        avg_ += val;
    }
    avg_ /= values_.size();
}

