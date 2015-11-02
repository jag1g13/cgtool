//
// Created by james on 25/08/15.
//

#include "GROInput.h"

#include <sstream>
#include <vector>
#include <sysexits.h>

#include <boost/algorithm/string.hpp>

using std::string;
using std::printf;
using std::vector;

/** \brief Contains data from a line of a GRO file. */
struct GROLine{
    int resnum = 0;
    std::string resname = "";
    std::string atomname = "";
    int atomnum = 0;
    float coords[3] = {0.f, 0.f, 0.f};
    float velocity[3] = {0.f, 0.f, 0.f};

    /** \brief Populate from string - line of GRO file. */
    int populate(const std::string &line){
        if(line.size() < 41) return 0;

        // Adjust for 7.2f formatting
        int float_len = 8;
        if(line.size() == 41) float_len = 7;

        resnum = std::stoi(line.substr(0, 5));
        resname = line.substr(5, 5);
        boost::trim(resname);
        atomname = line.substr(10, 5);
        boost::trim(atomname);
        atomnum = std::stoi(line.substr(15, 5));
        coords[0] = std::stof(line.substr(20, float_len));
        coords[1] = std::stof(line.substr(20 + float_len, float_len));
        coords[2] = std::stof(line.substr(20 + 2*float_len, float_len));

        if(line.size() >= 68){
            velocity[0] = std::stof(line.substr(44, float_len));
            velocity[1] = std::stof(line.substr(52, float_len));
            velocity[2] = std::stof(line.substr(60, float_len));
            return 2;
        }
        return 1;
    };

    /** \brief Print the parsed line for debugging. */
    void print() const{
        printf("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",
               resnum, resname.c_str(), atomname.c_str(),
               atomnum, coords[0], coords[1], coords[2]);
    }

    /** \brief Copy assignment operator. */
    GROLine& operator=(const GROLine &other){
        resnum = other.resnum;
        resname = other.resname;
        atomname = other.atomname;
        atomnum = other.atomnum;
        for(int i=0; i<3; i++){
            coords[i] = other.coords[i];
            velocity[i] = other.velocity[i];
        }

        return *this;
    }
};

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