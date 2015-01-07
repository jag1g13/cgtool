#include "cg_map.h"

#include <fstream>
#include <iostream>     // only needed for testing, get rid of this when done

#include <boost/algorithm/string.hpp>

#include "parser.h"

#define DEBUG true

using std::cout;
using std::endl;

CGMap::CGMap(){
}

CGMap::CGMap(string filename){
    fromFile(filename);
}

void CGMap::fromFile(string filename){
    string section;
    vector<string> substrs;
    Parser parser(filename);
    while(parser.getLine(&section, &substrs)){
        BeadMap new_bead;
        new_bead.cg_bead = substrs[0];
        new_bead.atoms = vector<string>(substrs.begin() + 1, substrs.end());
        new_bead.num_atoms = int(new_bead.atoms.size());
        mapping_.push_back(new_bead);
    }
    num_beads = mapping_.size();
    if(DEBUG){
        for(auto &bead : mapping_){
            std::cout << bead.cg_bead << " contains";
            for(auto &atom : bead.atoms){
                std::cout << " " << atom;
            }
            std::cout << std::endl;
        }
    }
}

void CGMap::initFrame(const Frame *aa_frame, Frame *cg_frame){
    cg_frame->num_ = aa_frame->num_;
    cg_frame->name_ = aa_frame->name_;
    cg_frame->prec_ = aa_frame->prec_;
    cg_frame->time_ = aa_frame->time_;
    //TODO finish copying values over
//    cg_frame->box_ = aa_frame->box_;
    int i = 0;
    for(auto &bead : mapping_) {
        cg_frame->atoms_.push_back(Atom(i));
        cg_frame->atoms_[i].atom_type_string = bead.cg_bead;
        cg_frame->atoms_[i].coords[0] = 0.f;
        cg_frame->atoms_[i].coords[1] = 0.f;
        cg_frame->atoms_[i].coords[2] = 0.f;
        cg_frame->name_to_num_.emplace(bead.cg_bead, i);
        cg_frame->num_to_name_.emplace(i, bead.cg_bead);
        for(auto &atomname : bead.atoms) {
            atomname_to_bead_.emplace(atomname, &bead);  // dictionary of atomnames to bead pointers
            //cout << bead->cg_bead << " contains " << *atomname << endl;
        }
//            for(auto &atom : aa_frame->atoms_){
        for(int j=0; j<aa_frame->numAtomsTrack_; j++) {
            for (string &atomname : bead.atoms){
//                cout << atomname << "\t" << aa_frame->atoms_[j].atom_type << endl;
                if(aa_frame->atoms_[j].atom_type == atomname){
                    cg_frame->atoms_[i].mass += aa_frame->atoms_[j].mass;
                    cg_frame->atoms_[i].charge += aa_frame->atoms_[j].charge;
//                    cout << "Match" << endl;
                }
            }
            if(aa_frame->atoms_[j].atom_type == bead.atoms[0]){
                bead.atom_nums.push_back(aa_frame->atoms_[j].atom_num);
            }
        }
        i++;
    }
    cg_frame->num_atoms_ = i+1;
    cg_frame->numAtomsTrack_ = i+1;

//    for(auto &atom : aa_frame->atoms_){
    for(int i=0; i<aa_frame->numAtomsTrack_; i++){
        // for atom in aa_frame that we care about
        const Atom *atom = &aa_frame->atoms_[i];
//        if(atom->atom_num < 20) cout << atomname_to_bead_[atom->atom_type] << endl;
        BeadMap *inbead = atomname_to_bead_[atom->atom_type];
//        if(atom->atom_num < 20) cout << &inbead << endl;
        inbead->mass += atom->mass;
        inbead->charge += atom->charge;
//        cout << inbead->mass << "\t" << inbead->charge << endl;
//        inbead->atom_nums.push_back(atom.atom_num);
//        this won't put them in the right order
    }

    cg_frame->isSetup_ = true;
    cout << "CG Frame" << endl;
    cg_frame->printAtoms();
    cout << "Done init cg_frame" << endl;
}

//TODO why is this a bool?
bool CGMap::apply(const Frame *aa_frame, Frame *cg_frame){
//    throw std::logic_error("Not implemented");
    bool status = true;
    if(!cg_frame->isSetup_) throw std::runtime_error("CG frame isn't setup");
    cg_frame->num_ = aa_frame->num_;
    cg_frame->time_ = aa_frame->time_;
//    cout << "Frame " << cg_frame->num_ << endl;
    if(mapType_ == MapType::ATOM){
        // if putting beads directly on the first atom in a bead
        int i = 0;
        for(auto &bead : mapping_){
//            int aa_num = aa_frame->name_to_num_[bead.atoms[0]];
            for(int j=0; j<3; j++){
                cg_frame->atoms_[i].coords[j] = aa_frame->atoms_[bead.atom_nums[0]].coords[j];
//                cout << cg_frame->atoms_[i].coords[j] << "\t";
            }
//            cout << endl;
            i++;
        }
    }else{
        throw std::runtime_error("Not implemented");
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

