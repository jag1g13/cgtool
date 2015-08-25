//
// Created by james on 25/08/15.
//

#include "GROInput.h"

#include <sstream>
#include <vector>
#include <sysexits.h>

using std::string;
using std::printf;
using std::vector;

GROInput::GROInput(const string &filename){
    if(openFile(filename)) throw std::runtime_error("Error opening GRO file");
}

GROInput::~GROInput(){
    closeFile();
}

int GROInput::openFile(const std::string &filename){
    file_.open(filename);

    string title;
    std::getline(file_, title);
    file_ >> natoms_;
    file_.seekg(0);

    if(!file_) return 1;
    return 0;
}

int GROInput::closeFile(){
    if(file_) file_.close();
    return 0;
}

int GROInput::readFrame(Frame &frame){
    // Discard top two lines of file - name and number of atoms
    string line;
    getline(file_, line);
    getline(file_, line);

    GROLine groline;
    for(int i=0; i<natoms_; i++){
        getline(file_, line);
        groline.populate(line);
//        if(i>0 && grolines[i].resname.compare(grolines[i-1].resname)) num_residues++;

        frame.atoms_[i].resnum = groline.resnum;
        frame.atoms_[i].atom_name = groline.atomname;
        frame.atoms_[i].coords[0] = groline.coords[0];
        frame.atoms_[i].coords[1] = groline.coords[1];
        frame.atoms_[i].coords[2] = groline.coords[2];
    }

    getline(file_, line);
    std::istringstream iss(line);
    iss >> frame.box_[0][0] >> frame.box_[1][1] >> frame.box_[2][2];

    return 0;
}

void GROInput::readResidues(vector<Residue> &residues){
    // Discard top two lines of file - name and number of atoms
    string line;
    getline(file_, line); getline(file_, line);

    // First pass to get number of residues
    int num_res = 0;
    GROLine current, prev;
    for(int i=0; i<natoms_; i++){
        std::getline(file_, line);
        current.populate(line);

        if(current.resname != prev.resname){
            num_res++;
        }
        prev = current;
    }

    if(residues.size() != num_res){
        printf("ERROR: Found %'d residue(s) in GRO and %'d in CFG\n",
               num_res, static_cast<int>(residues.size()));
        exit(EX_CONFIG);
    }

}