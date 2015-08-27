#include "frame.h"

#include <sstream>

#include <math.h>
#include <assert.h>
#include <sysexits.h>

#include "parser.h"
#include "small_functions.h"
#include "XTCOutput.h"
#include "GROOutput.h"
#include "LammpsDataOutput.h"
#include "LammpsTrjOutput.h"
#include "XTCInput.h"
#include "GROInput.h"

using std::string;
using std::vector;
using std::endl;
using std::stof;
using std::printf;
using std::map;

Frame::Frame(const Frame &frame, vector<Residue> *residues){
    num_ = frame.num_;
    name_ = frame.name_;
    time_ = frame.time_;
    step_ = frame.step_;
    boxType_ = frame.boxType_;

    if(residues == nullptr){
        residues_ = frame.residues_;
    }else{
        residues_ = residues;
    }

    isSetup_ = true;
}

Frame::Frame(const string &xtcname, const string &groname,
             vector<Residue> *residues){
    residues_ = residues;

    if(!initFromGRO(groname)){
        printf("ERROR: Something went wrong with reading GRO file\n");
        exit(EX_UNAVAILABLE);
    };

    trjIn_ = new XTCInput(xtcname);
}

Frame::~Frame(){
    isSetup_ = false;
    if(trjIn_) delete trjIn_;
    if(trjOut_) delete trjOut_;
}

void Frame::setupOutput(string xtcname, string topname){
    if(xtcname == "") xtcname = (*residues_)[0].resname + ".xtc";
    if(topname == "") topname = (*residues_)[0].resname + ".top";

    trjOut_ = new XTCOutput(numAtoms_, xtcname);
//    trjOut_ = new LammpsTrjOutput(numAtoms_, xtcname);
    writeTOP(topname);
    outputSetup_ = true;
}

void Frame::writeTOP(const string &filename){
    backup_old_file(filename);

    std::ofstream top(filename);
    if(!top.is_open()) throw std::runtime_error("Could not open output TOP file");

    top << "; Include forcefield parameters" << endl;
    top << "#include \"" << (*residues_)[0].resname << ".itp\"" << endl << endl;
    top << "[ system ]" << endl;
    top << (*residues_)[0].resname << endl << endl;
    top << "[ molecules ]" << endl;
    top << (*residues_)[0].resname << "\t\t" << (*residues_)[0].num_residues << endl;

    top.close();
}

bool Frame::outputTrajectoryFrame(){
    if(!outputSetup_) throw std::logic_error("Output has not been setup");

    return trjOut_->writeFrame(*this) == 0;
}

bool Frame::initFromGRO(const string &groname){
    GROInput in(groname);
    numAtoms_ = in.getNumAtoms();
    createAtoms(numAtoms_);
    in.readFrame(*this);
    in.readResidues(*residues_);

    // Create diciontary of residues
    for(Residue &res : *residues_){
        for(int i=0; i<res.num_atoms; i++){
            const int atom = res.start + i;
            res.name_to_num.insert(std::pair<string, int>(atoms_[atom].atom_name, i));
        }

        if(res.ref_atom_name != ""){
            if(res.name_to_num.find(res.ref_atom_name) == res.name_to_num.end()){
                printf("ERROR: Residue %6s does not contain reference atom %6s\n",
                       res.resname.c_str(), res.ref_atom_name.c_str());
                exit(EX_CONFIG);
            }else {
                res.ref_atom = res.name_to_num.at(res.ref_atom_name);
            }
        }
    }

    return true;
}

void Frame::createAtoms(int natoms){
    atoms_.resize(natoms);
    atomHas_.created = true;
}

void Frame::pbcAtom(int natoms){
    if(natoms < 0) natoms = numAtoms_;
    for(int i=0; i<natoms; i++){
        Atom &atom = atoms_[i];

        // For each coordinate wrap around into box
        for(int j=0; j<3; j++){
            while(atom.coords[j] < 0.) atom.coords[j] += box_[j][j];
            while(atom.coords[j] > box_[j][j]) atom.coords[j] -= box_[j][j];
        }
    }
}

void Frame::initFromITP(const string &itpname){
    // Require that atoms have been created
    assert(atomHas_.created);

    // Process topology file
    vector<string> substrs;
    Parser itp_parser(itpname, FileFormat::GROMACS);

    // How many atoms are there?  Per residue?  In total?
    if((*residues_)[0].num_atoms < 0){
        while(itp_parser.getLineFromSection("atoms", substrs, 4)){
            // Loop through all atoms in residue and take the last number
            if(substrs[3] == (*residues_)[0].resname) (*residues_)[0].num_atoms++;
        }
    }
    (*residues_)[0].calc_total();

    for(int i = 0; i < (*residues_)[0].num_atoms; i++){
        // Read data from topology file for each atom
        itp_parser.getLineFromSection("atoms", substrs, 5);
        const string type = substrs[1];
        const string name = substrs[4];
        atomHas_.atom_type = true;
        atomHas_.atom_name = true;

        double charge = 0.;
        if(substrs.size() >= 7){
            charge = atof(substrs[6].c_str());
            atomHas_.charge = true;
        }

        double mass = 1.;
        if(substrs.size() >= 8){
            mass = atof(substrs[7].c_str());
            atomHas_.mass = true;
        }

        for(int j = 0; j < (*residues_)[0].num_residues; j++){
            const int num = i + j * (*residues_)[0].num_atoms;
            atoms_[num] = Atom();
            atoms_[num].atom_type = type;
            atoms_[num].atom_name = name;
            atoms_[num].charge = charge;
            atoms_[num].mass = mass;
        }
    }

}

void Frame::initFromFLD(const std::string &fldname){
    // Require that atoms have been created and assigned names
    assert(atomHas_.created);
    assert(atomHas_.atom_type);

    Parser parser(fldname, FileFormat::GROMACS);

    map<string, double> c06;
    map<string, double> c12;

    vector<string> tokens;
    while(parser.getLineFromSection("atomtypes", tokens, 7)){
        c06[tokens[0]] = stof(tokens[5]);
        c12[tokens[0]] = stof(tokens[6]);
    }

    for(int i=0; i<(*residues_)[0].total_atoms; i++){
        atoms_[i].c06 = c06.at(atoms_[i].atom_type);
        atoms_[i].c12 = c12.at(atoms_[i].atom_type);
    }

    atomHas_.lj = true;
}

bool Frame::readNext(){
    return trjIn_->readFrame(*this) == 0;
}

void Frame::printAtoms(int natoms) const{
    assert(isSetup_);
    if(natoms == -1) natoms = numAtoms_;
    printf("  Num Name    Mass  Charge    Posx    Posy    Posz\n");
    for(int i=0; i<natoms; i++){
        const Atom &atom = atoms_[i];
        printf("%5i%5s%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",
               i, atom.atom_name.c_str(),
               atom.mass, atom.charge,
               atom.coords[0], atom.coords[1], atom.coords[2],
               atom.dipole[0], atom.dipole[1], atom.dipole[2], atom.dipole[3]);
    }
}

void Frame::outputSingleFrame(string filename) const{
    if(filename == "") filename = (*residues_)[0].resname + ".gro";
    GROOutput output(numAtoms_, filename);
//    LammpsDataOutput output(numAtoms_, filename);

    output.writeFrame(*this);
}

void Frame::printBox() const{
    // Print box vectors - assume cubic
    for(int i=0; i<3; i++){
        printf("%8.4f", box_[i][i]);
    }
    printf("\n");
}

