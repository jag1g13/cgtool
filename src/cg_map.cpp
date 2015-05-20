#include "cg_map.h"

#include <fstream>
#include <iostream>
#include <assert.h>

#include "parser.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;

CGMap::CGMap(const vector<Residue> &residues, const string &filename){
    aa_residue_ = residues[0];
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
        cout << "Could not find requested mapping type - assuming GC" << endl;
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
        new_bead.num_atoms = int(new_bead.atoms.size());

        mapping_.push_back(new_bead);

        // Print for debugging
        std::printf("%6s %6s %3i:", new_bead.name.c_str(), new_bead.type.c_str(), new_bead.num_atoms);
        for(const string &atom : new_bead.atoms) std::printf(" %s", atom.c_str());
        cout << endl;
    }
    numBeads_ = mapping_.size();
}

void CGMap::initFrame(const Frame &aa_frame, Frame &cg_frame){
    aa_residue_.print();

    // Create Frame and copy copyable data
    cg_residue_.resname = aa_residue_.resname;
    cg_residue_.start = 0;
    cg_residue_.num_atoms = numBeads_;
    cg_residue_.num_residues = aa_residue_.num_residues;
    cg_residue_.calc_total();
    cg_residue_.populated = true;

    cg_frame.residues_.resize(1);
    cg_frame.residues_[0] = cg_residue_;
    cg_frame.numAtoms_ = cg_residue_.total_atoms;
    cg_frame.atoms_.resize(aa_residue_.num_residues * aa_residue_.num_atoms);

    // Check if we have masses if CM mapping was requested
    if(mapType_ == MapType::CM && !aa_frame.atomHas_.mass){
        cout << "Centre of Mass mapping requires atom masses from ITP" << endl;
        cout << "Defaulting to Geometric Centre instead" << endl;
        mapType_ = MapType::GC;
    }

    // Create atom for each CG bead
    int i = 0;
    cg_frame.atomHas_.mass = true;
    for(BeadMap &bead : mapping_) {
        // Add bead to dictionaries so we can find it by name
        cg_frame.nameToNum_[bead.name] = i;

        // Calculate bead properties from atomistic frame
        for(const string &atomname : bead.atoms) {
            for(int j=aa_residue_.start; j<aa_residue_.start+aa_residue_.num_atoms; j++){
                if(aa_frame.atoms_[j].atom_name == atomname){
                    bead.mass += aa_frame.atoms_[j].mass;
                    bead.charge += aa_frame.atoms_[j].charge;
                    bead.c06 += aa_frame.atoms_[j].c06;
                    bead.c12 += aa_frame.atoms_[j].c12;
                    bead.atom_nums.push_back(j);
                }
            }
        }

        // Put properties into CG frame
        for(int j=0; j < aa_residue_.num_residues; j++){
            const int num_cg = i + j * cg_frame.residues_[0].num_atoms;
            cg_frame.atoms_[num_cg].atom_type = bead.name;
            cg_frame.atoms_[num_cg].atom_name = bead.name;
            cg_frame.atoms_[num_cg].charge = bead.charge;
            cg_frame.atoms_[num_cg].mass = bead.mass;
            cg_frame.atoms_[num_cg].resnum = j;
            cg_frame.atoms_[num_cg].c06 = bead.c06;
            cg_frame.atoms_[num_cg].c12 = bead.c12;
        }
        i++;

        // Check if
        if(bead.charge != 0.) cg_frame.atomHas_.charge = true;
        if(bead.mass == 0.) cg_frame.atomHas_.mass = false;
        if(bead.c06 != 0. && bead.c12 != 0.) cg_frame.atomHas_.lj = true;
    }

    cg_frame.numAtoms_ = i * aa_residue_.num_residues;

    apply(aa_frame, cg_frame);

    cg_frame.isSetup_ = true;
    cg_frame.atomHas_.atom_type = true;
    cg_frame.atomHas_.atom_name = true;
    cg_frame.atomHas_.resnum = true;
    cg_frame.atomHas_.coords = true;

    cout << "Done init cg_frame" << endl;
}

bool CGMap::apply(const Frame &aa_frame, Frame &cg_frame){
    bool status = true;
    if(!cg_frame.isSetup_) throw std::logic_error("CG frame isn't setup");
    cg_frame.num_ = aa_frame.num_;
    cg_frame.time_ = aa_frame.time_;
    cg_frame.step_ = aa_frame.step_;

    // Which mapping are we using?
    switch(mapType_){
        case MapType::ATOM:
            // If putting beads directly on the first atom in a bead
            for(int i = 0; i < mapping_.size(); i++){
                for(int j=0; j < aa_residue_.num_residues; j++){
                    const int num_cg = i + j*cg_frame.residues_[0].num_atoms;
                    const int num_aa = mapping_[i].atom_nums[0] + j*aa_frame.residues_[0].num_atoms;
                    cg_frame.atoms_[num_cg].coords[0] = aa_frame.atoms_[num_aa].coords[0];
                    cg_frame.atoms_[num_cg].coords[1] = aa_frame.atoms_[num_aa].coords[1];
                    cg_frame.atoms_[num_cg].coords[2] = aa_frame.atoms_[num_aa].coords[2];
                }
            }
            break;

        case MapType::GC:
            // Put bead at geometric centre of atoms
            for(int i = 0; i < mapping_.size(); i++){
                for(int j=0; j < aa_residue_.num_residues; j++){
                    const int num_cg = i + j*cg_frame.residues_[0].num_atoms;
                    cg_frame.atoms_[num_cg].coords[0] = 0.;
                    cg_frame.atoms_[num_cg].coords[1] = 0.;
                    cg_frame.atoms_[num_cg].coords[2] = 0.;

                    for(int k = 0; k < mapping_[i].num_atoms; k++){
                        const int num_aa = mapping_[i].atom_nums[k] + j*aa_frame.residues_[0].num_atoms;
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
            // Require that masses have been input
            assert(aa_frame.atomHas_.mass);

            // Put bead at centre of mass of atoms
            for(int i = 0; i < mapping_.size(); i++){
                for(int j = 0; j < aa_residue_.num_residues; j++){
                    const int num_cg = i + j * cg_frame.residues_[0].num_atoms;
                    cg_frame.atoms_[num_cg].coords[0] = 0.;
                    cg_frame.atoms_[num_cg].coords[1] = 0.;
                    cg_frame.atoms_[num_cg].coords[2] = 0.;

                    for(int k = 0; k < mapping_[i].num_atoms; k++){
                        const int num_aa = mapping_[i].atom_nums[k] + j*aa_frame.residues_[0].num_atoms;
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

