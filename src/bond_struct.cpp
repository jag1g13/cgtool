#include "bond_struct.h"

#include <stdexcept>
#include <array>

#include "small_functions.h"

using std::array;

BondStruct::BondStruct(const BondType type) : type_(type) {
    atomNums_.resize(static_cast<int>(type_));
}

double BondStruct::bondLength(const Frame &frame, const int offset) const{
    const int a = atomNums_[0] + offset;
    const int b = atomNums_[1] + offset;

    array<double, 3> vec;
    vec[0] = frame.atoms_[a].coords[0] - frame.atoms_[b].coords[0];
    vec[1] = frame.atoms_[a].coords[1] - frame.atoms_[b].coords[1];
    vec[2] = frame.atoms_[a].coords[2] - frame.atoms_[b].coords[2];

    return abs(vec);
}

double BondStruct::bondAngle(const Frame &frame, const int offset) const{
    const int a = atomNums_[0] + offset;
    const int b = atomNums_[1] + offset;
    int c, d;

    switch(type_){
        case BondType::DIHEDRAL:
            c = atomNums_[2] + offset;
            d = atomNums_[3] + offset;
            break;
        case BondType::ANGLE:
            c = atomNums_[1] + offset;
            d = atomNums_[2] + offset;
            break;
        default:
            throw std::logic_error("Passing a bond length as an angle");
    }

    array<double, 3> vec1, vec2;
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
