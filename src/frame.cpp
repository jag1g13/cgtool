#include "frame.h"

#include <iostream>

#include <math.h>
#include <assert.h>
#include <sstream>

#include "xdrfile_xtc.h"

#include "parser.h"

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
    atoms_.reserve(uint(natoms));
}

Frame::Frame(const Frame &frame){
    num_ = frame.num_;
    name_ = frame.name_;
    prec_ = frame.prec_;
    time_ = frame.time_;
    step_ = frame.step_;
    boxType_ = BoxType::CUBIC;
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            box_[i][j] = frame.box_[i][j];
            if(box_[i][j] > EPSILON) boxType_ = BoxType::TRICLINIC;
            // must be cubic box
//            if(i != j) assert(box_[i][j] < EPSILON);
        }
    }
}

Frame::Frame(const std::string topname, const std::string xtcname){
    setupFrame(topname, xtcname);
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

void Frame::setupOutput(const string &xtcname, const string &topname){
    char mode[2] = {'r', 'w'};
    if(xtcOutput_ == NULL) xtcOutput_ = xdrfile_open(xtcname.c_str(), &mode[1]);
    if(xtcOutput_ == NULL) throw std::runtime_error("Could not open XTC output");
    if(x_ == NULL){
        x_ = (rvec *) malloc(numAtoms_ * sizeof(*x_));
    }
    if(x_ == NULL) throw std::runtime_error("Couldn't allocate memory");

    std::ofstream top(topname);
    if(!top.is_open()) throw std::runtime_error("Could not open output TOP file");
    outputSetup_ = true;
}

bool Frame::writeToXtc(){
    if(!outputSetup_) throw std::runtime_error("Output has not been setup");
    // need to put atomic coordinates back into x_
    // either it's a CG frame and x_ is empty, or it's atomistic but may have been recentred
    for(int i=0; i<numAtoms_; i++){
        x_[i][0] = float(atoms_[i].coords[0]);
        x_[i][1] = float(atoms_[i].coords[1]);
        x_[i][2] = float(atoms_[i].coords[2]);
    }
    return exdrOK == write_xtc(xtcOutput_, numAtoms_, step_, time_, box_, x_, prec_);
}


bool Frame::setupFrame(const string &topname, const string &xtcname){
    if(isSetup_) throw std::logic_error("Frame has already been setup");
    char mode[2] = {'r', 'w'};
    int status = read_xtc_natoms(xtcname.c_str(), &numAtoms_);
    if(status != exdrOK){
        cout << "Could not open input XTC file" << endl;
        exit(-1);
    }
    xtcInput_ = xdrfile_open(xtcname.c_str(), &mode[0]);
//    if(output) xtc_out = xdrfile_open("out.xtc".c_str(), &mode[1]);
    num_ = 0;

    // init system from XTC file - libxdrfile library
    x_ = (rvec *)malloc(numAtoms_ * sizeof(*x_));
    status = read_xtc(xtcInput_, numAtoms_, &step_, &time_, box_, x_, &prec_);

//    recentreBox(0);

    // print box vectors
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
    while(top_parser.getLineFromSection("atoms", substrs)){
        numAtomsTrack_ = stoi(substrs[0]);
    }
    cout << numAtomsTrack_ << " atoms found in TOP/ITP" << endl;

    atoms_.resize(numAtomsTrack_);

    for(int i=0; i<numAtomsTrack_; i++){
        // read data from topology file for each atom
        // internal atom name is the res # and atom name from top/gro
        top_parser.getLineFromSection("atoms", substrs);
        string name = substrs[2] + substrs[4];

        atoms_[i] = Atom(i);
        atoms_[i].atom_type = name;
        atoms_[i].atom_num = int(stoi(substrs[0])-1);
        atoms_[i].resname = substrs[3];
        atoms_[i].charge = float(atof(substrs[6].c_str()));
        atoms_[i].mass = float(atof(substrs[7].c_str()));

        nameToNum_.emplace(atoms_[i].atom_type, i);

        atoms_[i].coords[0] = x_[i][0];
        atoms_[i].coords[1] = x_[i][1];
        atoms_[i].coords[2] = x_[i][2];
    }

    if(status == exdrOK) isSetup_ = true;
    printAtoms(numAtomsTrack_);
    return isSetup_;
}

bool Frame::readNext(){
    /**
    * \brief Read a frame from the XTC file into an existing Frame object
    *
    * Reads a frame into a pre-setup Frame object.
    * The same Frame object should be used for each frame to save time in allocation.
    */
    assert(isSetup_);
    invalid_ = false;
    int status = read_xtc(xtcInput_, numAtoms_, &step_, &time_, box_, x_, &prec_);
    if(status != exdrOK){
        xdrfile_close(xtcInput_);
        return status == exdrOK;
    }
//    recentreBox(0);
    for(int i = 0; i < numAtomsTrack_; i++){
        // overwrite coords of atoms stored in the current Frame
        atoms_[i].coords[0] = x_[i][0];
        atoms_[i].coords[1] = x_[i][1];
        atoms_[i].coords[2] = x_[i][2];
    }
    num_++;
    return status == exdrOK;
}

//TODO this doesn't solve the problem - we need all molecules to be whole
void Frame::recentreBox(const int atom_num){
    assert(isSetup_);
    assert(atom_num < numAtoms_);
//    assert(boxType_ == BoxType::CUBIC);

    float res_centre[3];

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

void Frame::printGRO(const string &filename, int natoms){
    assert(isSetup_);
    if(natoms == -1) natoms = numAtomsTrack_;
    FILE *gro = std::fopen(filename.c_str(), "w");
    if(gro == nullptr){
        cout << "Could not open gro file for writing" << endl;
        exit(-1);
    }
    std::stringstream stream;
    stream << "Generated by CGTOOL : " << name_ << "\n" << natoms << "\n";
    fprintf(gro, "%s", stream.str().c_str());
    double max[3], min[3];
    for(int i=0; i < natoms; i++){
        fprintf(gro, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
                1, "CG", atoms_[i].atom_type.c_str(), i+1,
                atoms_[i].coords[0], atoms_[i].coords[1], atoms_[i].coords[2]);
        for(int j=0; j < 3; j++){
            min[j] = fmin(min[j], atoms_[i].coords[j]);
            max[j] = fmax(max[j], atoms_[i].coords[j]);
        }
    }
    fprintf(gro, "%10.5f%10.5f%10.5f\n", max[0]-min[0], max[1]-min[1], max[2]-min[2]);
    fclose(gro);
}

double Frame::bondLength(const int a, const int b){
    return sqrt(pow((atoms_[a].coords[0] - atoms_[b].coords[0]), 2) +
            pow((atoms_[a].coords[1] - atoms_[b].coords[1]), 2) +
            pow((atoms_[a].coords[2] - atoms_[b].coords[2]), 2));
}

double Frame::bondLength(BondStruct &bond) {
    int a = bond.atomNums_[0];
    int b = bond.atomNums_[1];
    return bondLength(a, b);
}

double Frame::bondAngle(const int a, const int b, const int c, const int d){
    double vec1[3], vec2[3], mag1=0.f, mag2=0.f, dot = 0.f, angle=0.f;
    for(int i = 0; i < 3; i++){
        vec1[i] = atoms_[b].coords[i] - atoms_[a].coords[i];
        vec2[i] = atoms_[d].coords[i] - atoms_[c].coords[i];
        dot += vec1[i] * vec2[i];
    }
    mag1 = sqrt(pow(vec1[0], 2) + pow(vec1[1], 2) + pow(vec1[2], 2));
    mag2 = sqrt(pow(vec2[0], 2) + pow(vec2[1], 2) + pow(vec2[2], 2));
    angle = acos(dot / (mag1 * mag2));
    // if angle is NaN
    if(angle != angle){
//        printf("%s %s\n", numToName_[a].c_str(), numToName_[b].c_str());
        return 0.f;
    }
    return (180.f - (angle * 180.f / M_PI));
}

//TODO move this and bondLength into bond_struct.cpp
double Frame::bondAngle(BondStruct &bond){
    int a = bond.atomNums_[0];
    int b = bond.atomNums_[1];
    int c = bond.atomNums_[2];
    if(bond.atomNums_.size() == 4){
        int d = bond.atomNums_[3];
        return bondAngle(a, b, c, d);
    }else{
        return bondAngle(a, b, b, c);
    }
}
