#include "frame.h"

#include <fstream>
#include <iostream>
#include <cstring>
#include <limits>

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
    return (int)atoms_.size();
}

bool Frame::writeToXtc(t_fileio *xtc){
    return (bool)write_xtc(xtc, num_atoms_, step_, time_, box_, x_, prec_);
}

float Frame::bondLength(int a, int b){
    return (float)sqrt(pow((atoms_[a].coords[0] - atoms_[b].coords[0]), 2) +
            pow((atoms_[a].coords[1] - atoms_[b].coords[1]), 2) +
            pow((atoms_[a].coords[2] - atoms_[b].coords[2]), 2));
}

float Frame::bondLength(BondStruct *bond) {
    int a = name_to_num_[bond->atom_names[0]];
    int b = name_to_num_[bond->atom_names[1]];
    return bondLength(a, b);
}

float Frame::bondAngle(int a, int b, int c, int d){
    float vec1[3], vec2[3], mag1, mag2, dot = 0, angle;
    for(int i = 0; i < 3; i++){
        vec1[i] = atoms_[b].coords[i] - atoms_[a].coords[i];
        vec2[i] = atoms_[d].coords[i] - atoms_[c].coords[i];
        dot += vec1[i] * vec2[i];
    }
    mag1 = (float)sqrt(pow(vec1[0], 2) + pow(vec1[1], 2) + pow(vec1[2], 2));
    mag2 = (float)sqrt(pow(vec2[0], 2) + pow(vec2[1], 2) + pow(vec2[2], 2));
    angle = (float)acos(dot / (mag1 * mag2));
    return (180.f - (angle * 180.f / (float)M_PI));
}

float Frame::bondAngle(BondStruct *bond){
    int a = name_to_num_[bond->atom_names[0]];
    int b = name_to_num_[bond->atom_names[1]];
    int c = name_to_num_[bond->atom_names[2]];
    if(bond->atom_names.size() == 4){
        int d = name_to_num_[bond->atom_names[3]];
        return bondAngle(a, b, c, d);
    }else{
        return bondAngle(a, b, b, c);
    }
}

bool Frame::setupFrame(const char *groname, const char *topname, t_fileio *xtc){
    /**
    * \brief Create Frame, allocate atoms and read in data from start of XTC file
    *
    * GROMACS read_first_xtc() gets data from the XTC file about the system.
    * This function uses this data to create a Frame object to process this data
    */
    char line[40];
    int ok = 0, gro_num_atoms;
    num_ = 0;
    //float atom_charge, atom_mass;
    std::ifstream gro;
    gmx_bool bOK = 0;
    Atom *atom;
    ok = read_first_xtc(xtc, &num_atoms_, &step_, &time_, box_, &x_, &prec_, &bOK);
    gro.open(groname);
    if(gro.is_open()){
        gro.getline(line, 40);              // first line of gro is the run name
        cout << line << endl;
        gro >> gro_num_atoms;               // second line is the number of atoms
        if(gro_num_atoms != num_atoms_){
            cout << "XTC num atoms:" << num_atoms_ << endl;
            cout << "GRO num atoms:" << gro_num_atoms << endl;
            cout << "Number of atoms declared in XTC file "
                    "is not the same as declared in GRO file" << endl;
            throw std::runtime_error("Number of atoms does not match");
        }else{
            cout << "Found " << num_atoms_ << " atoms" << endl;
        }
        //allocateAtoms(num_atoms_);
        //Residue *res;
        string res_name_new, res_name_last;
        int res_loc = -1;
        for(int i = 0; i < num_atoms_; i++){       // now we can read the atoms
            atoms_.push_back(Atom(i));
            //atom = &(atoms_[i]);
            gro >> res_name_new >> atoms_[i].atom_type >> atoms_[i].atom_num;
            atoms_[i].atom_type[3] = '\0';
            gro.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            atoms_[i].atom_num--;
            atoms_[i].resid = res_name_new;
            if(res_name_new.compare(res_name_last) != 0){
                residues_.push_back(Residue(res_name_new));
                res_loc++;
                if(res_loc < 10 && res_loc > 0){
                    cout << "res: " << res_loc-1 << " resname: "<< residues_[res_loc-1].res_name;
                    cout << " size: " << residues_[res_loc-1].atoms.size() << endl;
                }
            }
            residues_[res_loc].atoms.push_back(atoms_[i].atom_num);
            residues_[res_loc].atom_names.push_back(atoms_[i].atom_type);
            memcpy(atoms_[i].coords, x_[i], 3 * sizeof(float));
            atoms_[i].charge = 0.f;
            atoms_[i].mass = 1.f;
            name_to_num_.emplace(atoms_[i].atom_type, i);
            num_to_name_.emplace(i, atoms_[i].atom_type);
            res_name_last = res_name_new;
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
    num_++;
    return ok && bOK;     //return True if it worked
}
