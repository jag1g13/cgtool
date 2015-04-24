#include "cg_map.h"

#include <fstream>
#include <iostream>

#include "parser.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;

CGMap::CGMap(const string &resname, const int numResidues, const string &filename){
    resname_ = resname;
    numResidues_ = numResidues;
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
        }else{
            cout << "Mapping type not recognised - assuming GC" << endl;
            mapType_ = MapType::GC;
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
    // Create Frame and copy copyable data
    cg_frame.numAtomsPerResidue_ = numBeads_;
    cg_frame.atoms_.resize(numResidues_ * cg_frame.numAtomsPerResidue_);

    // Create atom for each CG bead
    //TODO change this to support residues not from start
    int i = resBlockStart_;
    for(BeadMap &bead : mapping_) {
        for(int j=0; j < numResidues_; j++){
            const int num_cg = i + j * cg_frame.numAtomsPerResidue_;
            cg_frame.atoms_[num_cg].atom_type = bead.name;
            cg_frame.atoms_[num_cg].coords[0] = 0.;
            cg_frame.atoms_[num_cg].coords[1] = 0.;
            cg_frame.atoms_[num_cg].coords[2] = 0.;
        }

        // Add bead to dictionaries so we can find it by name
        cg_frame.nameToNum_.emplace(bead.name, i);
        for(const string &atomname : bead.atoms) {
            atomname_to_bead_.emplace(atomname, &bead);

            for(int j=0; j<aa_frame.numAtomsPerResidue_; j++){
                if(aa_frame.atoms_[j].atom_type == atomname){
                    cg_frame.atoms_[i].mass += aa_frame.atoms_[j].mass;
                    cg_frame.atoms_[i].charge += aa_frame.atoms_[j].charge;
                    bead.atom_nums.push_back(j);
                }
            }
        }

        // Copy values back into beads
        mapping_[i].mass = cg_frame.atoms_[i].mass;
        mapping_[i].charge = cg_frame.atoms_[i].charge;
        for(int j=0; j < numResidues_; j++){
            const int num_cg = i + j * cg_frame.numAtomsPerResidue_;
            cg_frame.atoms_[num_cg].mass = mapping_[i].mass;
            cg_frame.atoms_[num_cg].charge = mapping_[i].charge;
        }
        i++;
    }

    // Total number of atoms could include solvent later, but doesn't yet
    cg_frame.numAtoms_ = (i - resBlockStart_) * numResidues_;
    cg_frame.numAtomsTrack_ = cg_frame.numAtoms_;

    cg_frame.isSetup_ = true;
    apply(aa_frame, cg_frame);
    cout << "Done init cg_frame" << endl;
}

bool CGMap::apply(const Frame &aa_frame, Frame &cg_frame){
    bool status = true;
    if(!cg_frame.isSetup_) throw std::runtime_error("CG frame isn't setup");
    cg_frame.num_ = aa_frame.num_;
    cg_frame.time_ = aa_frame.time_;
    cg_frame.step_ = aa_frame.step_;

    // Remove 'invalid' marker - for frames where molecule crosses PBC
    cg_frame.invalid_ = false;

    // Which mapping are we using?
    switch(mapType_){
        case MapType::ATOM:
            // If putting beads directly on the first atom in a bead
            for(int i = 0; i < mapping_.size(); i++){
                for(int j=0; j < numResidues_; j++){
                    const int num_cg = i + j*cg_frame.numAtomsPerResidue_;
                    const int num_aa = mapping_[i].atom_nums[0] + j*aa_frame.numAtomsPerResidue_;
                    cg_frame.atoms_[num_cg].coords[0] = aa_frame.atoms_[num_aa].coords[0];
                    cg_frame.atoms_[num_cg].coords[1] = aa_frame.atoms_[num_aa].coords[1];
                    cg_frame.atoms_[num_cg].coords[2] = aa_frame.atoms_[num_aa].coords[2];
                }
            }
            break;

        case MapType::GC:
            // Put bead at geometric centre of atoms
            for(int i = 0; i < mapping_.size(); i++){
                for(int j=0; j < numResidues_; j++){
                    const int num_cg = i + j*cg_frame.numAtomsPerResidue_;
                    cg_frame.atoms_[num_cg].coords[0] = 0.;
                    cg_frame.atoms_[num_cg].coords[1] = 0.;
                    cg_frame.atoms_[num_cg].coords[2] = 0.;

                    for(int k = 0; k < mapping_[i].num_atoms; k++){
                        const int num_aa = mapping_[i].atom_nums[k] + j*aa_frame.numAtomsPerResidue_;
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
                for(int j = 0; j < numResidues_; j++){
                    const int num_cg = i + j * cg_frame.numAtomsPerResidue_;
                    cg_frame.atoms_[num_cg].coords[0] = 0.;
                    cg_frame.atoms_[num_cg].coords[1] = 0.;
                    cg_frame.atoms_[num_cg].coords[2] = 0.;

                    for(int k = 0; k < mapping_[i].num_atoms; k++){
                        const int num_aa = mapping_[i].atom_nums[k] + j*aa_frame.numAtomsPerResidue_;
                        const float mass = aa_frame.atoms_[num_aa].mass;
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

