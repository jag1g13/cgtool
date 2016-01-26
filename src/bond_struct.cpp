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
    const int c = atomNums_[2] + offset;

    array<double, 3> vec1, vec2;
//    vec1 = frame.atoms_[b].coords - frame.atoms_[a].coords;
    vec1[0] = frame.atoms_[b].coords[0] - frame.atoms_[a].coords[0];
    vec1[1] = frame.atoms_[b].coords[1] - frame.atoms_[a].coords[1];
    vec1[2] = frame.atoms_[b].coords[2] - frame.atoms_[a].coords[2];

//    vec2 = frame.atoms_[c].coords - frame.atoms_[b].coords;
    vec2[0] = frame.atoms_[c].coords[0] - frame.atoms_[b].coords[0];
    vec2[1] = frame.atoms_[c].coords[1] - frame.atoms_[b].coords[1];
    vec2[2] = frame.atoms_[c].coords[2] - frame.atoms_[b].coords[2];

    // Ensures angles between 0 and 360
    const double angle = M_PI - (acos(dot(vec1, vec2) / (abs(vec1) * abs(vec2))));
    return ((angle > 0 ? angle : (2*M_PI + angle)) * 180. / (M_PI));
}

double BondStruct::bondDihedral(const Frame &frame, const int offset) const{
    const int a = atomNums_[0] + offset;
    const int b = atomNums_[1] + offset;
    const int c = atomNums_[2] + offset;
    const int d = atomNums_[3] + offset;

    array<double, 3> vec1, vec2, vec3;
//    vec1 = frame.atoms_[b].coords - frame.atoms_[a].coords;
    vec1[0] = frame.atoms_[b].coords[0] - frame.atoms_[a].coords[0];
    vec1[1] = frame.atoms_[b].coords[1] - frame.atoms_[a].coords[1];
    vec1[2] = frame.atoms_[b].coords[2] - frame.atoms_[a].coords[2];

//    vec2 = frame.atoms_[c].coords - frame.atoms_[b].coords;
    vec2[0] = frame.atoms_[c].coords[0] - frame.atoms_[b].coords[0];
    vec2[1] = frame.atoms_[c].coords[1] - frame.atoms_[b].coords[1];
    vec2[2] = frame.atoms_[c].coords[2] - frame.atoms_[b].coords[2];

//    vec3 = frame.atoms_[d].coords - frame.atoms_[c].coords;
    vec3[0] = frame.atoms_[d].coords[0] - frame.atoms_[c].coords[0];
    vec3[1] = frame.atoms_[d].coords[1] - frame.atoms_[c].coords[1];
    vec3[2] = frame.atoms_[d].coords[2] - frame.atoms_[c].coords[2];

    array<double, 3> crossa, crossb;
    cross(vec1, vec2, crossa);
    cross(vec2, vec3, crossb);

    const double angle = M_PI - (acos(dot(crossa, crossb) / (abs(crossa) * abs(crossb))));
    return ((angle > 0 ? angle : (2*M_PI + angle)) * 180. / (M_PI));
}
