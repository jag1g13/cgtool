#include "frame.h"

#include <fstream>
#include <iostream>
#include <cstring>

#include <math.h>

using std::string;
using std::vector;
using std::cout;
using std::endl;

Frame::Frame(int num, int num_atoms, string name){
    name_ = name;
    step_ = num;
    num_atoms_ = num_atoms;
    atoms_.reserve(num_atoms);
}

Frame::Frame(const Frame* base_frame){
    name_ = base_frame->name_;
    step_ = base_frame->step_;
}

int Frame::allocateAtoms(int num_atoms){
    num_atoms_ = num_atoms;
    atoms_.reserve(num_atoms_);
    return 0;
}

bool Frame::writeToXtc(t_fileio *xtc){
    bool ok = write_xtc(xtc, num_atoms_, step_, time_, box_, x_, prec_);
    return ok;
    //throw std::logic_error("Not implemented");
}

float Frame::bondLength(int a, int b){
    return sqrt(pow((atoms_[a].coords[0] - atoms_[b].coords[0]), 2) +
            pow((atoms_[a].coords[1] - atoms_[b].coords[1]), 2) +
            pow((atoms_[a].coords[2] - atoms_[b].coords[2]), 2));
}

float Frame::bondAngle(int a, int b, int c, int d){
    float vec1[3], vec2[3], mag1, mag2, dot, angle;
    for(int i = 0; i < 3; i++){
        vec1[i] = atoms_[b].coords[i] - atoms_[a].coords[i];
        vec2[i] = atoms_[d].coords[i] - atoms_[c].coords[i];
        dot += vec1[i] * vec2[i];
    }
    mag1 = sqrt(pow(vec1[0], 2) + pow(vec1[1], 2) + pow(vec1[2], 2));
    mag2 = sqrt(pow(vec2[0], 2) + pow(vec2[1], 2) + pow(vec2[2], 2));
    angle = acos(dot / (mag1 * mag2));
    return 180. - (angle * 180. / M_PI);
}

bool Frame::setupFrame(char *groname, t_fileio *xtc){
    /**
    * \brief Create Frame, allocate atoms and read in data from start of XTC file
    *
    * GROMACS read_first_xtc() gets data from the XTC file about the system.
    * This function uses this data to create a Frame object to process this data
    */
    char res_name[5], line[40];
    int ok = 0, gro_num_atoms;
    //float atom_charge, atom_mass;
    std::ifstream gro;
    gmx_bool bOK = 0;
    Atom *atom;
    ok = read_first_xtc(xtc, &num_atoms_, &step_, &time_, box_, &x_, &prec_, &bOK);
    gro.open(groname);
    if(gro.is_open()){
        gro.getline(line, 40);              // first line of gro is the run name
        cout << line << endl;
        gro >> gro_num_atoms;                    // second line is the number of atoms
        if(gro_num_atoms != num_atoms_){
            cout << "XTC num atoms:" << num_atoms_ << endl;
            cout << "GRO num atoms:" << gro_num_atoms << endl;
            cout << "Number of atoms declared in XTC file "
                    "is not the same as declared in GRO file" << endl;
            throw std::runtime_error("Number of atoms does not match");
        }else{
            cout << "Found " << num_atoms_ << " atoms" << endl;
        }
        allocateAtoms(num_atoms_);
        for(int i = 0; i < num_atoms_; i++){       // now we can read the atoms
            atom = &(atoms_[i]);
            gro >> res_name >> atom->atom_type >> atom->atom_num;
            memcpy(atom->coords, x_[i], 3 * sizeof(float));
            atom->charge = 0.;
            atom->mass = 1.;
            name_to_num_.emplace(atom->atom_type, i);
            num_to_name_.emplace(i, atom->atom_type);
        }
        gro.close();
    } else{
        cout << "GRO file cannot be opened" << endl;
        throw std::runtime_error("Could not open GRO file");
    }
    return ok && bOK;                           // return True if it worked
}

bool Frame::readNext(t_fileio *xtc){
    /**
    * \brief Read a frame from the XTC file into an existing Frame object
    *
    * Reads a frame into a pre-setup Frame object.
    * The same Frame object should be used for each frame to save time in allocation.
    */
    int ok = 0, bOK = 0;
    //ok_out = write_xtc(xtc_out, *natoms, *step, *time, box, *x, *prec);
    ok = read_next_xtc(xtc, num_atoms_, &step_, &time_, box_, x_, &prec_, &bOK);
    for(int i = 0; i < num_atoms_; i++){
        memcpy(atoms_[i].coords, x_[i], 3 * sizeof(float)); // copy coordinates into an existing Atom
    }
    return ok && bOK;     //return True if it worked
}
