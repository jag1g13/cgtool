#include "bond_struct.h"

#include <iostream>

using std::cout;
using std::endl;

int BondStruct::totalNum_ = 0;

BondStruct::BondStruct(const int size){
    // keep track of serial number
    num_ = totalNum_;
    totalNum_++;
    atomNums_.resize(size);
    switch(size){
        case 2:
            type_ = BondType::LENGTH;
            break;
        case 3:
            type_ = BondType::ANGLE;
            break;
        case 4:
            type_ = BondType::DIHEDRAL;
            break;
    }
}

BondStruct::BondStruct(const BondStruct &other){
    unsigned long size = other.atomNums_.size();
    atomNums_.resize(size);
    for(unsigned int i = 0; i < size; i++){
        atomNums_[i] = other.atomNums_[i];
    }
}

void BondStruct::calcAvg(){
    for(const double &val : values_){
        avg_ += val;
    }
    avg_ /= values_.size();
}
