#include "cg_map.h"

#include <fstream>
#include <iostream>
#include <assert.h>
#include <cmath>

#include "parser.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;

CGMap::CGMap(const vector<Residue> *aa_residues, vector<Residue> *cg_residues,
             const string &filename){
    aaResidues_ = aa_residues;
    cgResidues_ = cg_residues;
    if(filename != "") fromFile(filename);
}

void CGMap::fromFile(const string &filename){
    // Which mapping type was requested - defaults to MapType::GC if not found
    vector<string> substrs;
    Parser parser(filename);
    if(parser.getLineFromSection("maptype", substrs, 1)){
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
    }else{
        printf("WARNING: Could not find requested mapping type - assuming GC\n");
        mapType_ = MapType::GC;
    }

    // Read in the bead mappings
    int i = 0;
    while(parser.getLineFromSection("mapping", substrs, 3)){
        BeadMap new_bead;
        new_bead.name = substrs[0];
        new_bead.type = substrs[1];
        new_bead.num = i++;
        new_bead.atoms = vector<string>(substrs.begin() + 2, substrs.end());
        new_bead.num_atoms = static_cast<int>(new_bead.atoms.size());

        mapping_.push_back(new_bead);

        // Print for debugging
        std::printf("%6s %6s %3i:", new_bead.name.c_str(), new_bead.type.c_str(), new_bead.num_atoms);
        for(const string &atom : new_bead.atoms) std::printf(" %s", atom.c_str());
        cout << endl;
    }
    numBeads_ = static_cast<int>(mapping_.size());
}

void CGMap::initFrame(const Frame &aa_frame, Frame &cg_frame){

    // Create Frame and copy copyable data
    (*cgResidues_).resize(1);
    (*cgResidues_)[0].resname = (*aaResidues_)[0].resname;
    (*cgResidues_)[0].start = 0;
    (*cgResidues_)[0].num_atoms = numBeads_;
    (*cgResidues_)[0].num_residues = (*aaResidues_)[0].num_residues;
    (*cgResidues_)[0].calc_total();
    (*cgResidues_)[0].populated = true;
    (*cgResidues_)[0].print();

    cg_frame.numAtoms_ = (*cgResidues_)[0].total_atoms;
    cg_frame.atoms_.resize(cg_frame.numAtoms_);

    // Check if we have masses if CM mapping was requested
    if(mapType_ == MapType::CM && !aa_frame.atomHas_.mass){
        printf("WARNING: Centre of Mass mapping requires atom masses from ITP\n");
        printf("WARNING: Defaulting to Geometric Centre instead\n");
        mapType_ = MapType::GC;
    }

    // Create atom for each CG bead
    int i = 0;
    cg_frame.atomHas_.mass = true;
    for(BeadMap &bead : mapping_) {
        // Add bead to dictionaries so we can find it by name
        (*cgResidues_)[0].name_to_num.insert(std::pair<string, int>(bead.name, i));

        // Vectors to store LJ values within a bead
        vector<int> c06s;
        vector<int> c12s;

        // Calculate bead properties from atomistic frame
        for(const string &atomname : bead.atoms){
            bool atom_found = false;
            for(int j=(*aaResidues_)[0].start; j<(*aaResidues_)[0].start+(*aaResidues_)[0].num_atoms; j++){
                if(aa_frame.atoms_[j].atom_name == atomname){
                    atom_found = true;
                    bead.mass += aa_frame.atoms_[j].mass;
                    bead.charge += aa_frame.atoms_[j].charge;
                    bead.atom_nums.push_back(j);

                    if(aa_frame.atomHas_.lj){
                        c06s.push_back(aa_frame.atoms_[j].c06);
                        c12s.push_back(aa_frame.atoms_[j].c12);
                    }
                }
            }

            if(!atom_found) printf("WARNING: Atom %s in bead %s not found\n",
                                   atomname.c_str(), bead.name.c_str());
        }

        if(aa_frame.atomHas_.lj){
            bead.c06 = calcLJ(c06s);
            bead.c12 = calcLJ(c12s);
            cg_frame.atomHas_.lj = true;
        }

        // Put properties into CG frame
        for(int j=0; j < (*aaResidues_)[0].num_residues; j++){
            const int num_cg = i + j * (*cgResidues_)[0].num_atoms;
            cg_frame.atoms_[num_cg].atom_type = bead.type;
            cg_frame.atoms_[num_cg].atom_name = bead.name;
            cg_frame.atoms_[num_cg].charge = bead.charge;
            cg_frame.atoms_[num_cg].mass = bead.mass;
            cg_frame.atoms_[num_cg].resnum = j;
            cg_frame.atoms_[num_cg].c06 = bead.c06;
            cg_frame.atoms_[num_cg].c12 = bead.c12;
        }
        i++;

        // Check if charge and mass are set for any atom
        if(bead.charge != 0.) cg_frame.atomHas_.charge = true;
        if(bead.mass == 0.) cg_frame.atomHas_.mass = false;
    }

    cg_frame.numAtoms_ = i * (*aaResidues_)[0].num_residues;

    apply(aa_frame, cg_frame);

    cg_frame.isSetup_ = true;
    cg_frame.atomHas_.atom_type = true;
    cg_frame.atomHas_.atom_name = true;
    cg_frame.atomHas_.resnum = true;
    cg_frame.atomHas_.coords = true;
}

double CGMap::calcLJ(const vector<int> &ljs){
    // Using GROMACS combination rule 1
    // A_b = sum_i( sum_j( sqrt(A_i*A_j) ) )

    double sum = 0.;
    for(const int a : ljs){
        for(const int b : ljs){
            sum += sqrt(a * b);
        }
    }

    return sum;
}

bool CGMap::apply(const Frame &aa_frame, Frame &cg_frame){
    bool status = true;
    if(!cg_frame.isSetup_) throw std::logic_error("CG frame isn't setup");
    cg_frame.num_ = aa_frame.num_;
    cg_frame.time_ = aa_frame.time_;
    cg_frame.step_ = aa_frame.step_;

    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            cg_frame.box_[i][j] = aa_frame.box_[i][j];
        }
    }


    // Which mapping are we using?
    switch(mapType_){
        case MapType::ATOM:
            // If putting beads directly on the first atom in a bead
            for(int i = 0; i < mapping_.size(); i++){
                for(int j=0; j < (*aaResidues_)[0].num_residues; j++){
                    const int num_cg = i + j*(*cgResidues_)[0].num_atoms;
                    const int num_aa = mapping_[i].atom_nums[0] + j*(*cgResidues_)[0].num_atoms;
                    cg_frame.atoms_[num_cg].coords[0] = aa_frame.atoms_[num_aa].coords[0];
                    cg_frame.atoms_[num_cg].coords[1] = aa_frame.atoms_[num_aa].coords[1];
                    cg_frame.atoms_[num_cg].coords[2] = aa_frame.atoms_[num_aa].coords[2];
                }
            }
            break;

        case MapType::GC:
            // Put bead at geometric centre of atoms
            for(int i = 0; i < mapping_.size(); i++){
                for(int j=0; j < (*aaResidues_)[0].num_residues; j++){
                    const int num_cg = i + j*(*cgResidues_)[0].num_atoms;
                    cg_frame.atoms_[num_cg].coords[0] = 0.;
                    cg_frame.atoms_[num_cg].coords[1] = 0.;
                    cg_frame.atoms_[num_cg].coords[2] = 0.;

                    for(int k = 0; k < mapping_[i].num_atoms; k++){
                        const int num_aa = mapping_[i].atom_nums[k] + j*(*aaResidues_)[0].num_atoms;
                        cg_frame.atoms_[num_cg].coords[0] += aa_frame.atoms_[num_aa].coords[0];
                        cg_frame.atoms_[num_cg].coords[1] += aa_frame.atoms_[num_aa].coords[1];
                        cg_frame.atoms_[num_cg].coords[2] += aa_frame.atoms_[num_aa].coords[2];
                    }
                    cg_frame.atoms_[num_cg].coords[0] /= mapping_[i].num_atoms;
                    cg_frame.atoms_[num_cg].coords[1] /= mapping_[i].num_atoms;
                    cg_frame.atoms_[num_cg].coords[2] /= mapping_[i].num_atoms;
                }
            }
            break;

        case MapType::CM:
            // Put bead at centre of mass of atoms
            for(int i = 0; i < mapping_.size(); i++){
                for(int j = 0; j < (*aaResidues_)[0].num_residues; j++){
                    const int num_cg = i + j * (*cgResidues_)[0].num_atoms;
                    cg_frame.atoms_[num_cg].coords[0] = 0.;
                    cg_frame.atoms_[num_cg].coords[1] = 0.;
                    cg_frame.atoms_[num_cg].coords[2] = 0.;

                    for(int k = 0; k < mapping_[i].num_atoms; k++){
                        const int num_aa = mapping_[i].atom_nums[k] + j*(*aaResidues_)[0].num_atoms;
                        const double mass = aa_frame.atoms_[num_aa].mass;
                        cg_frame.atoms_[num_cg].coords[0] += aa_frame.atoms_[num_aa].coords[0] * mass;
                        cg_frame.atoms_[num_cg].coords[1] += aa_frame.atoms_[num_aa].coords[1] * mass;
                        cg_frame.atoms_[num_cg].coords[2] += aa_frame.atoms_[num_aa].coords[2] * mass;
                    }
                    cg_frame.atoms_[num_cg].coords[0] /= mapping_[i].mass;
                    cg_frame.atoms_[num_cg].coords[1] /= mapping_[i].mass;
                    cg_frame.atoms_[num_cg].coords[2] /= mapping_[i].mass;
                }
            }
            break;
    }
    return status;
}

void CGMap::calcDipoles(const Frame &aa_frame, Frame &cg_frame){
    // For each molecule
    for(int k=0; k<(*cgResidues_)[0].num_residues; k++){
        // For each bead in the molecule
        for(int i = 0; i < numBeads_; i++){
            const BeadMap &bead_type = mapping_[i];
            Atom &cg_atom = cg_frame.atoms_[i];
            cg_atom.dipole[0] = cg_atom.dipole[1] = cg_atom.dipole[2] = 0.;

            // For each atom in the bead
            for(const int j : bead_type.atom_nums){
                const Atom &aa_atom = aa_frame.atoms_[j];
                // Rescale charges so bead charge is zero
                // This is how GMX_DIPOLE does it
                const double charge = aa_atom.charge -
                                      (cg_atom.charge * aa_atom.mass / bead_type.mass);
                cg_atom.dipole[0] += aa_atom.coords[0] * charge;
                cg_atom.dipole[1] += aa_atom.coords[1] * charge;
                cg_atom.dipole[2] += aa_atom.coords[2] * charge;
            }

            // Calculate magnitude
            cg_atom.dipole[3] = sqrt(cg_atom.dipole[0] * cg_atom.dipole[0] +
                                     cg_atom.dipole[1] * cg_atom.dipole[1] +
                                     cg_atom.dipole[2] * cg_atom.dipole[2]);
        }
    }
}
