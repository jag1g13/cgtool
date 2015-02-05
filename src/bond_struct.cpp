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

void BondStruct::boltzmannInversion(){
}

void BondStruct::binHistogram(const int bins){
    float max = float(avg_), min = float(avg_);
    histogram_.init(bins);

    for(const float val : values_){
        if(val < min) min = val;
        if(val > max) max = val;
    }

    float step = (max - min) / (bins-1);

    for(const float val : values_){
        int loc = int((val - min) / step);
        if(loc < 0 || loc > bins-1) cout << loc << endl;
        histogram_(loc)++;
    }
}