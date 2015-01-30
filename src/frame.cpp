#include "frame.h"

#include <iostream>

#include <math.h>
#include <assert.h>

#include "parser.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::stoi;
using std::printf;

Frame::Frame(int num, int natoms, string name){
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
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            box_[i][j] = frame.box_[i][j];
        }
    }
}

//Frame &Frame::operator=(const Frame &frame){
//    num_ = frame.num_;
//    name_ = frame.name_;
//    prec_ = frame.prec_;
//    time_ = frame.time_;
//    step_ = frame.step_;
//    for(int i=0; i<3; i++){
//        for(int j=0; j<3; j++){
//            box_[i][j] = frame.box_[i][j];
//        }
//    }
//    return *this;
//}

Frame::Frame(const std::string topname, const std::string xtcname){
    char mode[2] = {'r', 'w'};
    xtcInput_ = open_xtc(xtcname.c_str(), &mode[0]);
//    if(output) xtc_out = open_xtc("out.xtc", &mode[1]);
    setupFrame(topname, xtcInput_);
}

//TODO finish move constructor
//Frame::Frame(Frame&& frame){
//    isSetup_ = frame.isSetup_;
//    xtcInput_ = frame.xtcInput_;
//    num_ = frame.num_;
//    step_ = frame.step_;
//    numAtoms_ = frame.numAtoms_;
//    numAtomsTrack_ = frame.numAtomsTrack_;
//    atoms_ = frame.atoms_;
//    residues_ = frame.residues_;
//    time_ = frame.time_;
//    prec_ = frame.prec_;
////    box_ = frame.box_;
//    x_ = frame.x_;
//    name_ = frame.name_;
//    numToName_ = frame.numToName_;
//    nameToNum_ = frame.nameToNum_;
//}

Frame::~Frame(){
    assert(isSetup_);
    isSetup_ = false;
    if(x_) free(x_);
    // these lines actually result in less memory being freed
//    if(xtcInput_ != nullptr) close_xtc(xtcInput_);
//    xtcInput_ = nullptr;
//    if(xtcOutput_ != nullptr) close_xtc(xtcOutput_);
//    xtcOutput_ = nullptr;
    // don't seem to do anything
//    vector<Atom>().swap(atoms_);
//    vector<Residue>().swap(residues_);
}

int Frame::allocateAtoms(int num_atoms){
    numAtoms_ = num_atoms;
    atoms_.reserve(numAtoms_);
    return (int)atoms_.size();
}

void Frame::setupOutput(const string xtcname, const string topname){
    throw std::runtime_error("Not implemented");
    char mode[2] = {'r', 'w'};
    if(xtcOutput_ == NULL) xtcOutput_ = open_xtc(xtcname.c_str(), &mode[1]);
    if(x_ == NULL) x_ = (rvec*)malloc(numAtoms_ * sizeof(rvec));
    if(x_ == NULL) throw std::runtime_error("Couldn't allocate memory");

    std::ofstream top(topname);
    if(!top.is_open()) throw std::runtime_error("Could not open output TOP file");
}

bool Frame::writeToXtc(){
    throw std::runtime_error("Not implemented");
    if(!outputSetup_) throw std::runtime_error("Output has not been setup");
    // need to put atomic coordinates back into x_
    // either it's a CG frame and x_ is empty, or it's atomistic but may have been recentred
    for(int i=0; i<numAtoms_; i++){
        memcpy(&(x_[i][0]), atoms_[i].coords, 3);
    }
    return (bool)write_xtc(xtcOutput_, numAtoms_, step_, time_, box_, x_, prec_);
}


bool Frame::setupFrame(const std::string topname, t_fileio *xtc){
    if(isSetup_) throw std::runtime_error("Frame has already been setup");
    num_ = 0;
    gmx_bool bOK = 0;
    // init system from XTC file - GROMACS library
    int ok = read_first_xtc(xtc, &numAtoms_, &step_, &time_, box_, &x_, &prec_, &bOK);
    // recentre on first atom
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
    while(top_parser.getLineFromSection("atoms", &substrs)){
        numAtomsTrack_ = stoi(substrs[0]);
    }
    cout << numAtomsTrack_ << " atoms found in TOP/ITP" << endl;

    atoms_.resize(numAtomsTrack_);

    for(int i=0; i<numAtomsTrack_; i++){
        // read data from topology file for each atom
        // internal atom name is the res # and atom name from top/gro
        top_parser.getLineFromSection("atoms", &substrs);
        string name = substrs[2] + substrs[4];

        atoms_[i] = Atom(i);
        atoms_[i].atom_type = name;
        atoms_[i].atom_num = stoi(substrs[0])-1;
        atoms_[i].resname = substrs[3];
        atoms_[i].charge = float(atof(substrs[6].c_str()));
        atoms_[i].mass = float(atof(substrs[7].c_str()));

        nameToNum_.emplace(atoms_[i].atom_type, i);
        numToName_.emplace(i, atoms_[i].atom_type);

        atoms_[i].coords[0] = x_[i][0];
        atoms_[i].coords[1] = x_[i][1];
        atoms_[i].coords[2] = x_[i][2];
    }

    if(ok && bOK) isSetup_ = true;
//    printAtoms(numAtomsTrack_);
    return isSetup_;
}

//bool Frame::readNext(t_fileio *xtc){
bool Frame::readNext(){
    /**
    * \brief Read a frame from the XTC file into an existing Frame object
    *
    * Reads a frame into a pre-setup Frame object.
    * The same Frame object should be used for each frame to save time in allocation.
    */
    int ok = 0, bOK = 0;
    assert(isSetup_);
    invalid_ = false;
    ok = read_next_xtc(xtcInput_, numAtoms_, &step_, &time_, box_, x_, &prec_, &bOK);
//    recentreBox(0);
    for(int i = 0; i < numAtomsTrack_; i++){
        // overwrite coords of atoms stored in the current Frame
        atoms_[i].coords[0] = x_[i][0];
        atoms_[i].coords[1] = x_[i][1];
        atoms_[i].coords[2] = x_[i][2];
    }
    num_++;
    if(!ok || !bOK) close_xtc(xtcInput_);
    return ok && bOK;
}

//TODO check this works consistently
void Frame::recentreBox(const int atom_num){
    float box_centre[3];
    float res_centre[3];
    float offset[3];

    assert(atom_num < numAtoms_);
    res_centre[0] = x_[atom_num][0];
    res_centre[1] = x_[atom_num][1];
    res_centre[2] = x_[atom_num][2];
//    printf("res_centre: %8.4f%8.4f%8.4f\n", res_centre[0], res_centre[1], res_centre[2]);
    for(int i=0; i<numAtoms_; i++){
        for(int j=0; j<3; j++){
            x_[i][j] -= res_centre[j] - box_[2][0];
            if(x_[i][j] < -box_[2][0]) x_[i][j] += box_[0][0];
        }
    }
}

void Frame::printAtoms(int n){
    assert(isSetup_);
    if(n == -1) n = numAtomsTrack_;
    printf("  Name   Mass    Charge   Posx    Posy    Posz\n");
    for(int i=0; i<n; i++){
        printf("%6s %7.3f %7.3f %7.4f %7.4f %7.4f\n",
               atoms_[i].atom_type.c_str(), atoms_[i].mass, atoms_[i].charge,
               atoms_[i].coords[0], atoms_[i].coords[1], atoms_[i].coords[2]);
    }
}

float Frame::bondLength(int a, int b){
    return (float)sqrt(pow((atoms_[a].coords[0] - atoms_[b].coords[0]), 2) +
            pow((atoms_[a].coords[1] - atoms_[b].coords[1]), 2) +
            pow((atoms_[a].coords[2] - atoms_[b].coords[2]), 2));
}

float Frame::bondLength(BondStruct *bond) {
    int a = nameToNum_[bond->atom_names[0]];
    int b = nameToNum_[bond->atom_names[1]];
    return bondLength(a, b);
}

float Frame::bondAngle(int a, int b, int c, int d){
    float vec1[3], vec2[3], mag1=0.f, mag2=0.f, dot = 0.f, angle=0.f;
    for(int i = 0; i < 3; i++){
        vec1[i] = atoms_[b].coords[i] - atoms_[a].coords[i];
        vec2[i] = atoms_[d].coords[i] - atoms_[c].coords[i];
        dot += vec1[i] * vec2[i];
    }
    mag1 = (float)sqrt(pow(vec1[0], 2) + pow(vec1[1], 2) + pow(vec1[2], 2));
    mag2 = (float)sqrt(pow(vec2[0], 2) + pow(vec2[1], 2) + pow(vec2[2], 2));
    angle = (float)acos(dot / (mag1 * mag2));
    // if NaN
    if(angle != angle){
        printf("%d %d %d %d\n", a, b, c, d);
        return 0.f;
    }
    return (180.f - (angle * 180.f / (float)M_PI));
}

float Frame::bondAngle(BondStruct *bond){
    int a = nameToNum_[bond->atom_names[0]];
    int b = nameToNum_[bond->atom_names[1]];
    int c = nameToNum_[bond->atom_names[2]];
    if(bond->atom_names.size() == 4){
        int d = nameToNum_[bond->atom_names[3]];
        return bondAngle(a, b, c, d);
    }else{
        return bondAngle(a, b, b, c);
    }
}
