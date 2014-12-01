#include <iostream>
#include <vector>
#include <stdexcept>

#include <math.h>

#include "frame.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

Frame::Frame(int num, int num_atoms, string name){
    this->name = name;
    this->step = num;
    this->num_atoms = num_atoms;
    this->atoms.reserve(num_atoms);
    //cout << this->name << endl;
}

Frame::Frame(const Frame* base_frame){
    name = base_frame->name;
    step = base_frame->step;
}

int Frame::allocate_atoms(int num_atoms){
    this->num_atoms = num_atoms;
    this->atoms.reserve(num_atoms);
    return 0;
}

//TODO move as much as possible into the class
bool Frame::write_to_xtc(t_fileio *xtc){
    throw std::logic_error("Not implemented");
    gmx_bool bOK = 0;
    rvec *x;
    bool ok = 0;
    //bool ok = write_xtc(xtc, num_atoms, step, time, box, *x, prec);
    return ok && bOK;
}

float Frame::bond_length(int a, int b){
    return sqrt(pow((atoms[a].coords[0] - atoms[b].coords[0]), 2) +
            pow((atoms[a].coords[1] - atoms[b].coords[1]), 2) +
            pow((atoms[a].coords[2] - atoms[b].coords[2]), 2));
}

float Frame::bond_angle(int a, int b, int c, int d){
    float vec1[3], vec2[3], mag1, mag2, dot, angle;
    for(int i = 0; i < 3; i++){
        vec1[i] = atoms[b].coords[i] - atoms[a].coords[i];
        vec2[i] = atoms[d].coords[i] - atoms[c].coords[i];
        dot += vec1[i] * vec2[i];
    }
    mag1 = sqrt(pow(vec1[0], 2) + pow(vec1[1], 2) + pow(vec1[2], 2));
    mag2 = sqrt(pow(vec2[0], 2) + pow(vec2[1], 2) + pow(vec2[2], 2));
    angle = acos(dot / (mag1 * mag2));
    return 180. - (angle * 180. / M_PI);
}
