#include "bond_struct.h"

#include <iostream>
#include <stdexcept>

using std::cout;
using std::endl;

BondStruct::BondStruct(const int size){
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
        default:
            throw std::logic_error("BondStruct created with unexpected size");
    }
}

BondStruct::BondStruct(const BondStruct &other){
    unsigned long size = other.atomNums_.size();
    atomNums_.resize(size);
    for(unsigned long i = 0; i < size; i++){
        atomNums_[i] = other.atomNums_[i];
    }
    type_ = other.type_;
}
