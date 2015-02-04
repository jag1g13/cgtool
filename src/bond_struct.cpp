#include "bond_struct.h"

BondStruct::BondStruct(const int size){
    atomNames_.resize(size);
    atomNums_.resize(size);
}

void BondStruct::calcAvg(){
    for(float &val : values_){
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

    float step = (max - min) / bins;

    for(const float val : values_){
        int loc = int((val - min) / step);
        histogram_(loc)++;
    }
}