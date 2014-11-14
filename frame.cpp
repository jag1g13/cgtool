#include <iostream>
#include <vector>
#include <math.h>

#include "frame.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

Frame::Frame(int num, int num_atoms, string name){
    this->name = name;
    this->num = num;
    this->num_atoms = num_atoms;
    atoms.reserve(num_atoms);
    //atoms = (Atom*)malloc(num_atoms * sizeof(Atom));
    //if(atoms == NULL){
    //    cout << "Couldn't allocate memory for atoms" << endl;
    //}else{
    //    cout << "Allocated " << num_atoms << " atoms" << endl;
    //}
    cout << this->name << endl;
}

float Frame::bond_length(int a, int b){
    //distance between atoms i and j
    return sqrt(pow((atoms[a].coords[0] - atoms[b].coords[0]), 2) +
                pow((atoms[a].coords[1] - atoms[b].coords[1]), 2) +
                pow((atoms[a].coords[2] - atoms[b].coords[2]), 2));
}

float Frame::bond_angle(int a, int b, int c, int d){
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
