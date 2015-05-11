#include "frame.h"

#include <iostream>
#include <sstream>
#include <string>
#include <cstring>

#include <math.h>
#include <assert.h>
#include <sysexits.h>

#include "xdrfile_xtc.h"

#include "parser.h"
#include "small_functions.h"

#define EPSILON 0.00001

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::stoi;
using std::printf;

Residue& Residue::operator=(const Residue other){
    resname = other.resname;
    ref_atom = other.ref_atom;
    num_atoms = other.num_atoms;
    num_residues = other.num_residues;
    total_atoms = other.total_atoms;

    return *this;
}

Frame::Frame(const int num, const int natoms, const string name){
    assert(natoms >= 0);
    name_ = name;
    step_ = num;
    numAtoms_ = natoms;
    atoms_.reserve(natoms);
}

Frame::Frame(const Frame &frame){
    num_ = frame.num_;
    name_ = frame.name_;
    prec_ = frame.prec_;
    time_ = frame.time_;
    step_ = frame.step_;
    residue_ = frame.residue_;
    boxType_ = frame.boxType_;
    isSetup_ = true;
}

Frame::Frame(const string &itpname, const string &xtcname,
             const string &groname, const Residue &residue,
             const bool do_itp){
    residue_ = residue;
    initFromXTC(xtcname);
    if(do_itp) setupFrame(itpname, xtcname, groname);
}

Frame::~Frame(){
    isSetup_ = false;
    if(x_) free(x_);
    if(xtcOutput_){
        xdrfile_close(xtcOutput_);
    }
    if(xtcInput_){
        xdrfile_close(xtcInput_);
    }
}

void Frame::setupOutput(string xtcname, string topname){
    char mode[2] = {'r', 'w'};
    if(xtcname == "") xtcname = residue_.resname + ".xtc";
    if(topname == "") topname = residue_.resname + ".top";
    backup_old_file(xtcname);
    backup_old_file(topname);
    if(!xtcOutput_) xtcOutput_ = xdrfile_open(xtcname.c_str(), &mode[1]);
    if(!xtcOutput_) throw std::runtime_error("Could not open XTC output");

    if(!x_) x_ = new rvec[numAtoms_];
    if(!x_) throw std::runtime_error("Couldn't allocate memory for XTC output");

    std::ofstream top(topname);
    if(!top.is_open()) throw std::runtime_error("Could not open output TOP file");
    top << "; Include forcefield parameters" << endl;
    top << "#include \"" << residue_.resname << "CG.itp\"" << endl << endl;
    top << "[ system ]" << endl;
    top << residue_.resname << endl << endl;
    top << "[ molecules ]" << endl;
    top << residue_.resname << "\t\t" << residue_.num_residues << endl;
    top.close();
    outputSetup_ = true;
}

bool Frame::writeToXtc(){
    if(!outputSetup_) throw std::logic_error("Output has not been setup");
    // need to put atomic coordinates back into x_
    // either it's a CG frame and x_ is empty, or it's atomistic but may have been recentred
    for(int i=0; i<numAtoms_; i++){
        x_[i][0] = float(atoms_[i].coords[0]);
        x_[i][1] = float(atoms_[i].coords[1]);
        x_[i][2] = float(atoms_[i].coords[2]);
    }
    return exdrOK == write_xtc(xtcOutput_, numAtoms_, step_, time_, box_, x_, prec_);
}

void Frame::openOtherXTC(const Frame &frame){
    if(xtcInput_) xdrfile_close(xtcInput_);
    xtcInput_ = frame.xtcInput_;
}

bool Frame::initFromXTC(const string &xtcname){
    // How many atoms?  Prepare Frame for reading
    int status = read_xtc_natoms(xtcname.c_str(), &numAtoms_);
    if(status != exdrOK){
        cout << "Could not open input XTC file" << endl;
        exit(EX_NOINPUT);
    }
    printf("XTC contains %'d atoms\n", numAtoms_);
    xtcInput_ = xdrfile_open(xtcname.c_str(), "r");
    num_ = 0;

    // Init system from XTC file - libxdrfile library
    x_ = (rvec *)malloc(numAtoms_ * sizeof(*x_));
    status = read_xtc(xtcInput_, numAtoms_, &step_, &time_, box_, x_, &prec_);

    // Print box vectors
    cout << "Box vectors" << endl;
    boxType_ = BoxType::CUBIC;
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            box_[i][j] = box_[i][j];
            if(i!=j && box_[i][j] > EPSILON) boxType_ = BoxType::TRICLINIC;
            // Require cubic box
//            if(i != j) assert(box_[i][j] < EPSILON);
        }
    }
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            printf("%8.4f", box_[i][j]);
        }
        cout << endl;
    }

    if(status == exdrOK) isSetup_ = true;
    return isSetup_;
}

void Frame::initFromGRO(const string &groname, vector<Residue> &residues){
    cout << "Using " << groname << endl;
    FILE *gro = fopen(groname.c_str(), "r");
    if(gro == NULL){
        printf("Could not open GRO file for reading\n");
        exit(EX_NOINPUT);
    }

    // Have we been given a complete list of residues or do we need to create it?
    bool populate_residues = false;
    if(residues.size() == 0) populate_residues = true;

    char system_name[80];
    fgets(system_name, 80, gro);
    printf("%s\n", system_name);

    int num_atoms = 0;
    fscanf(gro, "%d", &num_atoms);
    printf("GRO contains %'d atoms\n", num_atoms);

    // The columns of a GRO file
    int resnum;
    char resname[6];
    char atomname[6];
    int atomnum;
    float atom_coords[3];
    float atom_vel[3];

    while(fscanf(gro, "%5d%5s%5s%5d%*8f%*8f%*8f%*8f%*8f%*8f",
                 &resnum, resname, atomname, &atomnum)){
        atoms_[atomnum-1] = Atom();
        atoms_[atomnum-1].atom_type = atomname;
        atoms_[atomnum-1].resnum = resnum;
        if(populate_residues){
            // Is it a new residue
            if(residues.back().resname != resname){

            }
        }
    }
    fclose(gro);
}

void Frame::copyCoordsIntoAtoms(){
    assert(atoms_.size() == numAtoms_);
    for(int i=0; i<numAtoms_; i++){
        atoms_[i].coords[0] = x_[i][0];
        atoms_[i].coords[1] = x_[i][1];
        atoms_[i].coords[2] = x_[i][2];
    }
}

void Frame::setupFrame(const string &itpname, const string &xtcname,
                       const string &groname){

    // Process topology file
    vector<string> substrs;
    Parser itp_parser(itpname, FileFormat::GROMACS);
    // How many atoms are there?  Per residue?  In total?
    while(itp_parser.getLineFromSection("atoms", substrs, 4)){
        // If we don't have a resname from cfg, be backward compatible
        // Loop through all atoms and take the last number
        if(residue_.resname != ""){
            if(substrs[3] == residue_.resname) residue_.num_atoms = stoi(substrs[0]);
        } else{
            residue_.num_atoms = stoi(substrs[0]);
        }
    }

    cout << "Found " << residue_.num_atoms << " atoms in ITP per " << residue_.resname << endl;
    numAtomsTrack_ = residue_.num_residues * residue_.num_atoms;
    atoms_.resize(numAtomsTrack_);

    for(int i = 0; i < residue_.num_atoms; i++){
        // read data from topology file for each atom
        // internal atom name is the res # and atom name from top/gro
        itp_parser.getLineFromSection("atoms", substrs, 5);
        const string name = substrs[4];
        double charge = 0.;
        if(substrs.size() >= 7) charge = atof(substrs[6].c_str());
        double mass = 1.;
        if(substrs.size() >= 8) mass = atof(substrs[7].c_str());
        //TODO why doesn't this work
        nameToNum_[atoms_[i].atom_type] = i;

        for(int j = 0; j < residue_.num_residues; j++){
            const int num = i + j * residue_.num_atoms;
            atoms_[num] = Atom();
            atoms_[num].atom_type = name;
            atoms_[num].charge = charge;
            atoms_[num].mass = mass;

            atoms_[num].coords[0] = x_[num][0];
            atoms_[num].coords[1] = x_[num][1];
            atoms_[num].coords[2] = x_[num][2];
        }
    }

}

bool Frame::readNext(){
    assert(isSetup_);
    int status = read_xtc(xtcInput_, numAtoms_, &step_, &time_, box_, x_, &prec_);
    if(status != exdrOK) return false;
    for(int i = 0; i < numAtomsTrack_; i++){
        // Overwrite coords of atoms stored in the current Frame
        atoms_[i].coords[0] = x_[i][0];
        atoms_[i].coords[1] = x_[i][1];
        atoms_[i].coords[2] = x_[i][2];
    }
    num_++;
    return true;
}

void Frame::recentreBox(const int atom_num){
    throw std::logic_error("Function not implemented");
    assert(isSetup_);
    assert(atom_num < numAtoms_);
    assert(boxType_ == BoxType::CUBIC);

    double res_centre[3];

    res_centre[0] = x_[atom_num][0];
    res_centre[1] = x_[atom_num][1];
    res_centre[2] = x_[atom_num][2];
    printf("res_centre: %8.4f%8.4f%8.4f\n", res_centre[0], res_centre[1], res_centre[2]);
//    for(int i=0; i<numAtoms_; i++){
//        for(int j=0; j<3; j++){
//            x_[i][j] -= res_centre[j];
//            if(x_[i][j] < -box_[j][j]) x_[i][j] += box_[j][j];
//        }
//    }
}


void Frame::printAtoms(int natoms){
    assert(isSetup_);
    if(natoms == -1) natoms = numAtomsTrack_;
    printf("   Num   Name   Mass   Charge   Posx    Posy    Posz\n");
    for(int i=0; i<natoms; i++){
        printf("%6i %6s %7.3f %7.3f %7.4f %7.4f %7.4f\n",
               nameToNum_.at(atoms_[i].atom_type), atoms_[i].atom_type.c_str(),
               atoms_[i].mass, atoms_[i].charge,
               atoms_[i].coords[0], atoms_[i].coords[1], atoms_[i].coords[2]);
    }
}

void Frame::printGRO(string filename, int natoms){
    assert(isSetup_);
    if(filename == "") filename = residue_.resname + ".gro";
    if(natoms == -1) natoms = numAtomsTrack_;
    backup_old_file(filename);

    FILE *gro = std::fopen(filename.c_str(), "w");
    if(gro == nullptr){
        cout << "Could not open gro file for writing" << endl;
        exit(EX_CANTCREAT);
    }

    // Set initial guess for system bounds
    double min[3], max[3];
    for(int j=0; j < 3; j++){
        min[j] = atoms_[0].coords[j];
        max[j] = atoms_[0].coords[j];
    }

    // Print data to GRO file
    std::stringstream stream;
    stream << "Generated by CGTOOL : " << residue_.resname << "\n" << natoms << "\n";
    fprintf(gro, "%s", stream.str().c_str());
    for(int i=0; i < natoms; i++){
        fprintf(gro, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
                1+(i/residue_.num_atoms), residue_.resname.c_str(), atoms_[i].atom_type.c_str(), i+1,
                atoms_[i].coords[0], atoms_[i].coords[1], atoms_[i].coords[2]);

        for(int j=0; j < 3; j++){
            min[j] = fmin(min[j], atoms_[i].coords[j]);
            max[j] = fmax(max[j], atoms_[i].coords[j]);
        }
    }
    fprintf(gro, "%10.5f%10.5f%10.5f\n", max[0]-min[0], max[1]-min[1], max[2]-min[2]);

    fclose(gro);
}

double Frame::bondLength(BondStruct &bond, const int offset) {
    const int a = bond.atomNums_[0] + offset;
    const int b = bond.atomNums_[1] + offset;

    double vec[3];
    vec[0] = atoms_[a].coords[0] - atoms_[b].coords[0];
    vec[1] = atoms_[a].coords[1] - atoms_[b].coords[1];
    vec[2] = atoms_[a].coords[2] - atoms_[b].coords[2];

    return abs(vec);
}

//TODO move this and bondLength into bond_struct.cpp
double Frame::bondAngle(BondStruct &bond, const int offset){
    const int a = bond.atomNums_[0] + offset;
    const int b = bond.atomNums_[1] + offset;
    int c, d;

    switch(bond.atomNums_.size()){
        case 4:
            c = bond.atomNums_[2] + offset;
            d = bond.atomNums_[3] + offset;
            break;
        case 3:
            c = bond.atomNums_[1] + offset;
            d = bond.atomNums_[2] + offset;
            break;
        default:
            throw std::logic_error("Passing a bond length as an angle");
    }

    double vec1[3], vec2[3];
    vec1[0] = atoms_[b].coords[0] - atoms_[a].coords[0];
    vec1[1] = atoms_[b].coords[1] - atoms_[a].coords[1];
    vec1[2] = atoms_[b].coords[2] - atoms_[a].coords[2];

    vec2[0] = atoms_[d].coords[0] - atoms_[c].coords[0];
    vec2[1] = atoms_[d].coords[1] - atoms_[c].coords[1];
    vec2[2] = atoms_[d].coords[2] - atoms_[c].coords[2];

    const double angle = acos(dot(vec1, vec2) / (abs(vec1) * abs(vec2)));
    return (180. - (angle * 180. / M_PI));
}
