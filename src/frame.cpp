#include "frame.h"

#include <iostream>
#include <sstream>
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
    resname_ = frame.resname_;
    numResidues_ = frame.numResidues_;
    boxType_ = BoxType::CUBIC;
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            box_[i][j] = frame.box_[i][j];
            if(i!=j && box_[i][j] > EPSILON) boxType_ = BoxType::TRICLINIC;
            // must be cubic box
//            if(i != j) assert(box_[i][j] < EPSILON);
        }
    }
}

Frame::Frame(const string &topname, const string &xtcname, const string &groname, const string &resname, const int numResidues){
    resname_ = resname;
    numResidues_ = numResidues;
    setupFrame(topname, xtcname, groname);
}

Frame::~Frame(){
    assert(isSetup_);
    isSetup_ = false;
    if(x_ != nullptr) free(x_);
    if(xtcOutput_ != nullptr){
        xdrfile_close(xtcOutput_);
        xtcOutput_ = NULL;
    }
    if(xtcInput_ != nullptr){
//        xdrfile_close(xtcInput_);
//        xtcInput_ = NULL;
    }
}

void Frame::setupOutput(string xtcname, string topname){
    char mode[2] = {'r', 'w'};
    if(xtcname == "") xtcname = resname_ + "CG.xtc";
    if(topname == "") topname = resname_ + "CG.top";
    backup_old_file(xtcname);
    backup_old_file(topname);
    if(!xtcOutput_) xtcOutput_ = xdrfile_open(xtcname.c_str(), &mode[1]);
    if(!xtcOutput_) throw std::runtime_error("Could not open XTC output");

    if(!x_) x_ = (rvec *) malloc(numAtoms_ * sizeof(*x_));
    if(!x_) throw std::runtime_error("Couldn't allocate memory for XTC output");

    std::ofstream top(topname);
    if(!top.is_open()) throw std::runtime_error("Could not open output TOP file");
    top << "; Include forcefield parameters" << endl;
    top << "#include \"" << resname_ << "CG.itp\"" << endl << endl;
    top << "[ system ]" << endl;
    top << resname_ << endl << endl;
    top << "[ molecules ]" << endl;
    top << resname_ << "\t\t" << numResidues_ << endl;
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

bool Frame::setupFrame(const string &topname, const string &xtcname, const string &groname){
    if(isSetup_) throw std::logic_error("Frame has already been setup");

    // How many atoms?  Prepare Frame for reading
    int status = read_xtc_natoms(xtcname.c_str(), &numAtoms_);
    if(status != exdrOK){
        cout << "Could not open input XTC file" << endl;
        exit(EX_NOINPUT);
    }
    xtcInput_ = xdrfile_open(xtcname.c_str(), "r");
    num_ = 0;

    // Init system from XTC file - libxdrfile library
    x_ = (rvec *)malloc(numAtoms_ * sizeof(*x_));
    status = read_xtc(xtcInput_, numAtoms_, &step_, &time_, box_, x_, &prec_);

    // Print box vectors
    cout << "Box vectors" << endl;
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            printf("%8.4f", box_[i][j]);
        }
        cout << endl;
    }

    // Process topology file
    vector<string> substrs;
    Parser top_parser(topname, ParserFormat::GROMACS);
    // How many atoms are there?  Per residue?  In total?
    while(top_parser.getLineFromSection("atoms", substrs)){
        // If we don't have a resname from cfg, be backward compatible
        // Loop through all atoms and take the last number
        if(resname_ != ""){
            if(substrs[3] == resname_) numAtomsPerResidue_ = stoi(substrs[0]);
        }else{
            numAtomsPerResidue_ = stoi(substrs[0]);
        }
    }

    cout << "Found " << numAtomsPerResidue_ << " atoms in ITP per " << resname_ << endl;
    numAtomsTrack_ = numResidues_ * numAtomsPerResidue_;
    atoms_.resize(numAtomsTrack_);

    // If we have a GRO then find the resname we're looking for
    FILE *gro = fopen(groname.c_str(), "r");
    char name[6];
    int start = 0;
    if(gro != NULL){
        while(fscanf(gro, "%*5d%5s%*5s%*5d%*8f%*8f%*8f%*8f%*8f%*8f", name)){
            if(strcmp(name, resname_.c_str()) == 0) break;
            start++;
        }
    }


    for(int i=0; i<numAtomsPerResidue_; i++){
        // read data from topology file for each atom
        // internal atom name is the res # and atom name from top/gro
        top_parser.getLineFromSection("atoms", substrs);
        string name = substrs[2] + substrs[4];
        const double charge = atof(substrs[6].c_str());
        const double mass = atof(substrs[7].c_str());
        nameToNum_.emplace(atoms_[i].atom_type, i);

        for(int j=0; j<numResidues_; j++){
            const int num = i + j*numAtomsPerResidue_;
            atoms_[num] = Atom();
            atoms_[num].atom_num = num;
            atoms_[num].atom_type = name;
            atoms_[num].charge = charge;
            atoms_[num].mass = mass;

            atoms_[num].coords[0] = x_[num+start][0];
            atoms_[num].coords[1] = x_[num+start][1];
            atoms_[num].coords[2] = x_[num+start][2];
        }
    }

    if(status == exdrOK) isSetup_ = true;
    return isSetup_;
}

bool Frame::readNext(){
    assert(isSetup_);
    invalid_ = false;
    int status = read_xtc(xtcInput_, numAtoms_, &step_, &time_, box_, x_, &prec_);
    if(status != exdrOK){
        xdrfile_close(xtcInput_);
        return status == exdrOK;
    }
    for(int i = 0; i < numAtomsTrack_; i++){
        // Overwrite coords of atoms stored in the current Frame
        atoms_[i].coords[0] = x_[i][0];
        atoms_[i].coords[1] = x_[i][1];
        atoms_[i].coords[2] = x_[i][2];
    }
    num_++;
    return status == exdrOK;
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
    for(int i=0; i<numAtoms_; i++){
        for(int j=0; j<3; j++){
//            x_[i][j] -= res_centre[j];
//            if(x_[i][j] < -box_[j][j]) x_[i][j] += box_[j][j];
        }
    }
}


void Frame::printAtoms(int natoms){
    assert(isSetup_);
    if(natoms == -1) natoms = numAtomsTrack_;
    printf("   Num   Name   Mass   Charge   Posx    Posy    Posz\n");
    for(int i=0; i<natoms; i++){
        printf("%6i %6s %7.3f %7.3f %7.4f %7.4f %7.4f\n",
               nameToNum_[atoms_[i].atom_type], atoms_[i].atom_type.c_str(),
               atoms_[i].mass, atoms_[i].charge,
               atoms_[i].coords[0], atoms_[i].coords[1], atoms_[i].coords[2]);
    }
}

void Frame::printGRO(string filename, int natoms){
    assert(isSetup_);
    if(filename == "") filename = resname_ + "CG.gro";
    if(natoms == -1) natoms = numAtomsTrack_;
    backup_old_file(filename);

    FILE *gro = std::fopen(filename.c_str(), "w");
    if(gro == nullptr){
        cout << "Could not open gro file for writing" << endl;
        exit(EX_CANTCREAT);
    }

    std::stringstream stream;
    stream << "Generated by CGTOOL : " << resname_ << "\n" << natoms << "\n";
    fprintf(gro, "%s", stream.str().c_str());
    double min[3]={1000., 1000., 1000.}, max[3] = {-1000., -1000., -1000.};
    for(int i=0; i < natoms; i++){
        fprintf(gro, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
                1+(i/numAtomsPerResidue_), resname_.c_str(), atoms_[i].atom_type.c_str(), i+1,
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
    for(int i = 0; i < 3; i++){
        vec1[i] = atoms_[b].coords[i] - atoms_[a].coords[i];
        vec2[i] = atoms_[d].coords[i] - atoms_[c].coords[i];
    }

    const double angle = acos(dot(vec1, vec2) / (abs(vec1) * abs(vec2)));
    return (180. - (angle * 180. / M_PI));
}
