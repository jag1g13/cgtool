#include "bond_struct.h"

#include <stdexcept>

#include "small_functions.h"

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
    size_t size = other.atomNums_.size();
    atomNums_.resize(size);
    for(size_t i = 0; i < size; i++){
        atomNums_[i] = other.atomNums_[i];
    }
    type_ = other.type_;
}

double BondStruct::bondLength(const Frame &frame, const int offset) {
    const int a = atomNums_[0] + offset;
    const int b = atomNums_[1] + offset;

    double vec[3];
    vec[0] = frame.atoms_[a].coords[0] - frame.atoms_[b].coords[0];
    vec[1] = frame.atoms_[a].coords[1] - frame.atoms_[b].coords[1];
    vec[2] = frame.atoms_[a].coords[2] - frame.atoms_[b].coords[2];

    return abs(vec);
}

double BondStruct::bondAngle(const Frame &frame, const int offset){
    const int a = atomNums_[0] + offset;
    const int b = atomNums_[1] + offset;
    int c, d;

    switch(atomNums_.size()){
        case 4:
            c = atomNums_[2] + offset;
            d = atomNums_[3] + offset;
            break;
        case 3:
            c = atomNums_[1] + offset;
            d = atomNums_[2] + offset;
            break;
        default:
            throw std::logic_error("Passing a bond length as an angle");
    }

    double vec1[3], vec2[3];
    vec1[0] = frame.atoms_[b].coords[0] - frame.atoms_[a].coords[0];
    vec1[1] = frame.atoms_[b].coords[1] - frame.atoms_[a].coords[1];
    vec1[2] = frame.atoms_[b].coords[2] - frame.atoms_[a].coords[2];

    vec2[0] = frame.atoms_[d].coords[0] - frame.atoms_[c].coords[0];
    vec2[1] = frame.atoms_[d].coords[1] - frame.atoms_[c].coords[1];
    vec2[2] = frame.atoms_[d].coords[2] - frame.atoms_[c].coords[2];

    // Ensures angles between 0 and 360
    const double angle = M_PI - (acos(dot(vec1, vec2) / (abs(vec1) * abs(vec2))));
    return ((angle > 0 ? angle : (2*M_PI + angle)) * 180. / (M_PI));
}
