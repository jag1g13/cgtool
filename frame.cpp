#include <iostream>
#include <vector>
#include <math.h>

#include "frame.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

Frame::Frame(int num, int num_atoms, string name){
    /**
    * \brief Create Frame passing frame number, number of atoms to store and the frame name
    *
    * If we don't know the number of atoms at creation
    * this can be set later using Frame::allocate_atoms()
    */
    this->name = name;
    this->num = num;
    this->num_atoms = num_atoms;
    this->atoms.reserve(num_atoms);
    //cout << this->name << endl;
}

int Frame::allocate_atoms(int num_atoms){
    /**
    * \brief Allocate space for a number of atoms
    *
    * Used if the number of atoms isn't known at time of creation
    */
    this->num_atoms = num_atoms;
    this->atoms.reserve(num_atoms);
    return 0;
}

float Frame::bond_length(int a, int b){
    /**
    * \brief Calculate distance between two atoms
    */
    return sqrt(pow((atoms[a].coords[0] - atoms[b].coords[0]), 2) +
                pow((atoms[a].coords[1] - atoms[b].coords[1]), 2) +
                pow((atoms[a].coords[2] - atoms[b].coords[2]), 2));
}

float Frame::bond_angle(int a, int b, int c, int d){
    /**
    * \brief Calculate angle between vectors a->b and c->d
    *
    * To be used for bond angles (b=c) and dihedrals (b=/=c)
    */
    float vec1[3], vec2[3], mag1, mag2, dot, angle;
    for(int i=0; i<3; i++){
        vec1[i] = atoms[b].coords[i] - atoms[a].coords[i];
        vec2[i] = atoms[d].coords[i] - atoms[c].coords[i];
        dot += vec1[i] * vec2[i];
    }
    mag1 = sqrt(pow(vec1[0], 2) + pow(vec1[1], 2) + pow(vec1[2], 2));
    mag2 = sqrt(pow(vec2[0], 2) + pow(vec2[1], 2) + pow(vec2[2], 2));
    angle = acos(dot / (mag1 * mag2));
    return 180. - (angle * 180. / M_PI);
}
