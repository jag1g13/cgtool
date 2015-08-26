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
    file_.seekg(0);
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

    frame.atomHas_.atom_name = true;
    frame.atomHas_.resnum = true;
    frame.atomHas_.coords = true;

    return 0;
}

void GROInput::readResidues(vector<Residue> &residues){
    // Discard top two lines of file - name and number of atoms
    file_.seekg(0);
    string line;
    getline(file_, line); getline(file_, line);

    // First pass to get number of residues
    int num_res = 0;
    GROLine current, prev;
    Residue *res = &(residues[0]);
    for(int i=0; i<natoms_; i++){
        std::getline(file_, line);
        current.populate(line);

        if(current.resname != prev.resname){
            num_res++;

            if(residues.size() < num_res){
                printf("ERROR: Not all residues are listed in CFG\n");
                exit(EX_CONFIG);
            }

            res = &(residues[num_res-1]);
            res->start = i;
            res->resname = current.resname;
            res->total_atoms = 0;
            res->num_residues = 0;

        }

        if(current.resnum != prev.resnum){
            res->num_residues++;
        }

        res->total_atoms++;
        prev = current;
    }


    if(residues.size() != num_res){
        printf("ERROR: Found %'d residue(s) in GRO and %'d in CFG\n",
               num_res, static_cast<int>(residues.size()));
        exit(EX_CONFIG);
    }

    // Second pass for missing values
    res = &(residues[num_res-1]);
    res->end = natoms_;
    for(int i=0; i<num_res-1; i++){
        res = &(residues[i]);
        Residue *res_next = &(residues[i+1]);

        res->end = res_next->start;
        res->set_num_atoms(res->total_atoms / res->num_residues);
        res->populated = true;
    }
    res = &(residues[num_res-1]);
    res->set_num_atoms(res->total_atoms / res->num_residues);
    res->populated = true;
}