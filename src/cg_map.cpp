#include "cg_map.h"

#include <fstream>
#include <iostream>     // only needed for testing, get rid of this when done

#include <boost/algorithm/string.hpp>

#include "parser.h"

#define DEBUG false

using std::cout;
using std::endl;

CGMap::CGMap(){
}

CGMap::CGMap(string filename){
    fromFile(filename);
}

void CGMap::fromFile(string filename){
    vector<string> substrs;
    Parser parser(filename);
    // which mapping type was requested - defaults to MapType::ATOM if not found
    while(parser.getLineFromSection("maptype", &substrs)){
        if(substrs[0] == "CM"){
            mapType_ = MapType::CM;
            cout << "Using CM mapping" << endl;
        }else if(substrs[0] == "GC"){
            mapType_ = MapType::GC;
            cout << "Using GC mapping" << endl;
        }else if(substrs[0] == "ATOM"){
            mapType_ = MapType::ATOM;
            cout << "Using ATOM mapping" << endl;
        }
    }
    // read in the bead mappings
    while(parser.getLineFromSection("mapping", &substrs)){
        BeadMap new_bead;
        new_bead.cg_bead = substrs[0];
        new_bead.atoms = vector<string>(substrs.begin() + 1, substrs.end());
        new_bead.num_atoms = int(new_bead.atoms.size());

        // print for debugging
        cout << new_bead.cg_bead << " " << new_bead.num_atoms << " ";
        for(auto atom : new_bead.atoms){
            cout << atom << " ";
        }
        cout << endl;
        mapping_.push_back(new_bead);
    }
    num_beads = int(mapping_.size());
}

Frame CGMap::initFrame(const Frame &aa_frame){
    Frame cg_frame(aa_frame);

    int i = 0;
    for(auto &bead : mapping_) {
        cg_frame.atoms_.push_back(Atom(i));
        cg_frame.atoms_[i].atom_type = bead.cg_bead;
        cg_frame.atoms_[i].coords[0] = 0.f;
        cg_frame.atoms_[i].coords[1] = 0.f;
        cg_frame.atoms_[i].coords[2] = 0.f;
        cg_frame.nameToNum_.emplace(bead.cg_bead, i);
        cg_frame.numToName_.emplace(i, bead.cg_bead);
        for(auto &atomname : bead.atoms) {
            // dictionary of atom to bead they're in
            atomname_to_bead_.emplace(atomname, &bead);
        }
        for(int j=0; j<aa_frame.numAtomsTrack_; j++) {
            for (string &atomname : bead.atoms){
                if(aa_frame.atoms_[j].atom_type == atomname){
                    cg_frame.atoms_[i].mass += aa_frame.atoms_[j].mass;
                    cg_frame.atoms_[i].charge += aa_frame.atoms_[j].charge;
                }
                if(aa_frame.atoms_[j].atom_type == atomname){
                    bead.atom_nums.push_back(aa_frame.atoms_[j].atom_num);
                }
            }
        }
        i++;
    }
    cg_frame.numAtoms_ = i;
    cg_frame.numAtomsTrack_ = i;

    for(int i=0; i<aa_frame.numAtomsTrack_; i++){
//        cout << "Mapping bead " << i << endl;
        // for atom in aa_frame that we care about
        const Atom *atom = &(aa_frame.atoms_[i]);
        if(atomname_to_bead_.count(atom->atom_type)){
            BeadMap *inbead = atomname_to_bead_[atom->atom_type];
            inbead->mass += atom->mass;
            inbead->charge += atom->charge;
        }else{
            cout << "Ignoring atom " << i << endl;
        }
    }

    cg_frame.isSetup_ = true;
    cg_frame.printAtoms();
    apply(aa_frame, cg_frame);
    cout << "CG Frame" << endl;
    cg_frame.printAtoms();
    cout << "Done init cg_frame" << endl;

    return cg_frame;
}

bool CGMap::apply(const Frame &aa_frame, Frame &cg_frame){
    bool status = true;
    if(!cg_frame.isSetup_) throw std::runtime_error("CG frame isn't setup");
    cg_frame.num_ = aa_frame.num_;
    cg_frame.time_ = aa_frame.time_;
    cg_frame.invalid_ = false;

    // which mapping are we using?
    switch(mapType_){
        case MapType::ATOM:
            // if putting beads directly on the first atom in a bead
            for(int i = 0; i < mapping_.size(); i++){
                for(int j = 0; j < 3; j++){
                    cg_frame.atoms_[i].coords[j] =
                            aa_frame.atoms_[mapping_[i].atom_nums[0]].coords[j];
                }
            }
            break;

        case MapType::GC:
            // put bead at geometric centre of atoms
            for(int i = 0; i < mapping_.size(); i++){
                cout << "Applying bead " << i << " Atom ";
                for(int j = 0; j < mapping_[i].num_atoms; j++){
                    int num = mapping_[i].atom_nums[j];
                    cout << num << " ";
                    cg_frame.atoms_[i].coords[0] += aa_frame.atoms_[num].coords[0];
                    cg_frame.atoms_[i].coords[1] += aa_frame.atoms_[num].coords[1];
                    cg_frame.atoms_[i].coords[2] += aa_frame.atoms_[num].coords[2];
                }
                cout << endl;
                for(int k=0; k < 3; k++){
                    cg_frame.atoms_[i].coords[k] /= mapping_[i].num_atoms;
                }
            }
            break;

        case MapType::CM:
            // put bead at centre of mass of atoms
            for(int i = 0; i < mapping_.size(); i++){
                for(int j = 0; j < mapping_[i].num_atoms; j++){
                    for(int k = 0; k < 3; k++){
                        cg_frame.atoms_[i].coords[k] +=
                                aa_frame.atoms_[mapping_[i].atom_nums[j]].coords[k] *
                                aa_frame.atoms_[mapping_[i].atom_nums[j]].mass;
                    }
                }
                for(int k=0; k < 3; k++){
                    cg_frame.atoms_[i].coords[k] /= mapping_[i].mass;
                }
            }
            break;
    }
    return status;
}

