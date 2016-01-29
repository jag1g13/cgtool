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

    array<double, 3> vec = frame.atoms_[b].coords - frame.atoms_[a].coords;
    pbcWrap(vec, frame.boxDiag_);

    return abs(vec);
}

double BondStruct::bondAngle(const Frame &frame, const int offset) const{
    const int a = atomNums_[0] + offset;
    const int b = atomNums_[1] + offset;
    const int c = atomNums_[2] + offset;

    array<double, 3> vec1 = frame.atoms_[b].coords - frame.atoms_[a].coords;
    array<double, 3> vec2 = frame.atoms_[c].coords - frame.atoms_[b].coords;
    pbcWrap(vec1, frame.boxDiag_);
    pbcWrap(vec2, frame.boxDiag_);

    return (M_PI - angle(vec1, vec2)) * 180. / M_PI;
}

double BondStruct::bondDihedral(const Frame &frame, const int offset) const{
    const int a = atomNums_[0] + offset;
    const int b = atomNums_[1] + offset;
    const int c = atomNums_[2] + offset;
    const int d = atomNums_[3] + offset;

    array<double, 3> vec1 = frame.atoms_[b].coords - frame.atoms_[a].coords;
    array<double, 3> vec2 = frame.atoms_[c].coords - frame.atoms_[b].coords;
    array<double, 3> vec3 = frame.atoms_[d].coords - frame.atoms_[c].coords;
    pbcWrap(vec1, frame.boxDiag_);
    pbcWrap(vec2, frame.boxDiag_);
    pbcWrap(vec3, frame.boxDiag_);

    array<double, 3> crossa, crossb, crossc;
    cross(vec1, vec2, crossa);
    cross(vec2, vec3, crossb);
    cross(crossa, crossb, crossc);

    double ang = angle(crossa, crossb) * 180. / M_PI;
    const double dir = dot(vec2, crossc);
    return dir < 0 ? ang : -ang;
}
