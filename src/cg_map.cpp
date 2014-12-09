#include "cg_map.h"

#include <fstream>
#include <iostream>     // only needed for testing, get rid of this when done

//#include <boost/tokenizer.hpp>
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
        mapping_.push_back(new_bead);
    }
    num_beads = mapping_.size();
    if(DEBUG){
        for(auto &i : mapping_){
            std::cout << i.cg_bead << " contains";
            for(auto &j : i.atoms){
                std::cout << " " << j;
            }
            std::cout << std::endl;
        }
    }
}

void CGMap::initFrame(Frame *aa_frame, Frame *cg_frame){
    for(auto &bead : mapping_){
        for(auto &atomname : bead){
            atomname_to_bead_.emplace(*atomname, bead);  // dictionary of atomnames to bead pointers
            //cout << bead->cg_bead << " contains " << *atomname << endl;
        }
    }
    for(auto &atom : aa_frame->atoms_){
        // for atom in aa_frame
        BeadMap* inbead = atomname_to_bead_[atom.atom_type];
        // need to work out bonding
    }
}

bool CGMap::apply(const Frame *aa_frame, Frame *cg_frame){
    throw std::logic_error("Not implemented");

    bool status = true;
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

