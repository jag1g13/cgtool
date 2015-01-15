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

    }
    // read in the bead mappings
    while(parser.getLineFromSection("mapping", &substrs)){
        BeadMap new_bead;
        new_bead.cg_bead = substrs[0];
        new_bead.atoms = vector<string>(substrs.begin() + 1, substrs.end());
        new_bead.num_atoms = int(new_bead.atoms.size());
        mapping_.push_back(new_bead);
    }
    num_beads = mapping_.size();
    if(DEBUG){
        for(auto &bead : mapping_){
//            std::cout << bead.cg_bead << " contains";
            printf("%s\n", bead.cg_bead.c_str());
            for(auto &atom : bead.atoms){
//                std::cout << " " << atom;
                printf("%6s", atom.c_str());
            }
            printf("\n");
//            std::cout << std::endl;
        }
    }
}

//void CGMap::initFrame(const Frame *aa_frame, Frame *cg_frame){
Frame CGMap::initFrame(const Frame &aa_frame){
//    cg_frame->name_ = aa_frame->name_;
//    cg_frame->prec_ = aa_frame->prec_;
//    cg_frame->num_ = aa_frame->num_;
//    cg_frame->time_ = aa_frame->time_;
//    cg_frame->step_ = aa_frame->step_;
//    for(int i=0; i<3; i++){
//        for(int j=0; j<3; j++){
//            cg_frame->box_[i][j] = aa_frame->box_[i][j];
//        }
//    }
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
            atomname_to_bead_.emplace(atomname, &bead);  // dictionary of atomnames to bead pointers
            //cout << bead->cg_bead << " contains " << *atomname << endl;
        }
//            for(auto &atom : aa_frame->atoms_){
        for(int j=0; j<aa_frame.numAtomsTrack_; j++) {
            for (string &atomname : bead.atoms){
//                cout << atomname << "\t" << aa_frame->atoms_[j].atom_type << endl;
                if(aa_frame.atoms_[j].atom_type == atomname){
                    cg_frame.atoms_[i].mass += aa_frame.atoms_[j].mass;
                    cg_frame.atoms_[i].charge += aa_frame.atoms_[j].charge;
//                    cout << "Match" << endl;
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

//    for(auto &atom : aa_frame->atoms_){
    for(int i=0; i<aa_frame.numAtomsTrack_; i++){
//        cout << "starting atom " << i << endl;
        // for atom in aa_frame that we care about
        const Atom *atom = &(aa_frame.atoms_[i]);
//        cout << atom->atom_type << "\t";
//        if(atom->atom_num < 20) cout << atomname_to_bead_[atom->atom_type] << endl;
        BeadMap *inbead = atomname_to_bead_[atom->atom_type];
//        if(atom->atom_num < 20) cout << &inbead << endl;
//        cout << inbead->cg_bead << endl;
        inbead->mass += atom->mass;
        inbead->charge += atom->charge;
//        cout << inbead->mass << "\t" << inbead->charge << endl;
//        inbead->atom_nums.push_back(atom.atom_num);
//        this won't put them in the right order
//        cout << "finishing atom" << endl;
    }

    cg_frame.isSetup_ = true;
    apply(aa_frame, cg_frame);
    cout << "CG Frame" << endl;
    cg_frame.printAtoms();
    cout << "Done init cg_frame" << endl;

    return cg_frame;
}

bool CGMap::apply(const Frame &aa_frame, Frame &cg_frame){
//    throw std::logic_error("Not implemented");
    bool status = true;
    if(!cg_frame.isSetup_) throw std::runtime_error("CG frame isn't setup");
    cg_frame.num_ = aa_frame.num_;
    cg_frame.time_ = aa_frame.time_;
//    cout << "Frame " << cg_frame->num_ << endl;
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
                for(int j = 0; j < mapping_[i].num_atoms; j++){
                    for(int k = 0; k < 3; k++){
                        cg_frame.atoms_[i].coords[k] +=
                                aa_frame.atoms_[mapping_[i].atom_nums[j]].coords[k];
                    }
                }
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

// the mapping function used in traj_process (Python)
/*def map_cg_solvent_within_loop(curr_frame, frame, cg_frame=0):
    """
    perform CG mapping using cg_map list of lists
            with current cg_map does a simple heavy atom mapping

    will be CM or GC depending on how 'Atom.mass' was set previously
    if mapping is changed in cg_map to include other atoms

            should remove the setup code into its own function (or the main xtc setup)
    """
    global cg_atom_nums
    if curr_frame == 0:
        cg_frame = Frame(curr_frame, cg_atom_nums)
    cg_frame.num = curr_frame
    for i, site in enumerate(cg_sites):
        coords = np.zeros(3)
        tot_mass = 0.
        charge = 0.
        for atom in cg_map[i]:
            mass = frame.atoms[sugar_atom_nums[atom]].mass
            tot_mass = tot_mass + mass
            coords = coords + mass*frame.atoms[sugar_atom_nums[atom]].loc
            charge = charge + frame.atoms[sugar_atom_nums[atom]].charge
        coords /= tot_mass  # number of atoms cancels out
        if curr_frame == 0:
            cg_frame.atoms.append(Atom(site, coords, charge))
        else:
            cg_frame.atoms[i] = Atom(site, coords, charge)
        if curr_frame == 0:
            cg_atom_nums[site] = i
    j = len(cg_sites)
    for atom in frame.atoms:
        if atom.atom_type == "OW":
            if curr_frame == 0:
                cg_frame.atoms.append(Atom("OW", atom.loc, 0.0))
            else:
                cg_frame.atoms[j] = Atom("OW", atom.loc, 0.0)
            j += 1
    return cg_frame*/

