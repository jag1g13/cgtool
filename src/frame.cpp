#include "frame.h"

#include <fstream>
#include <iostream>
#include <cstring>
#include <limits>
#include <exception>

#include <math.h>
#include <assert.h>

#include "parser.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;

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

Frame::Frame(const std::string groname, const std::string topname, const std::string cfgname, const std::string xtcname){
    char mode[2] = {'r', 'w'};
    xtcInput_ = open_xtc(xtcname.c_str(), &mode[0]);
//    if(output) xtc_out = open_xtc("out.xtc", &mode[1]);
    setupFrame(groname, topname, cfgname, xtcInput_);
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
    float vec1[3], vec2[3], mag1, mag2, dot = 0, angle;
    for(int i = 0; i < 3; i++){
        vec1[i] = atoms_[b].coords[i] - atoms_[a].coords[i];
        vec2[i] = atoms_[d].coords[i] - atoms_[c].coords[i];
        dot += vec1[i] * vec2[i];
    }
    mag1 = (float)sqrt(pow(vec1[0], 2) + pow(vec1[1], 2) + pow(vec1[2], 2));
    mag2 = (float)sqrt(pow(vec2[0], 2) + pow(vec2[1], 2) + pow(vec2[2], 2));
    angle = (float)acos(dot / (mag1 * mag2));
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

bool Frame::setupFrame(const std::string groname, const std::string topname, const std::string cfgname, t_fileio *xtc){
    if(isSetup_) throw std::runtime_error("Frame has already been setup");
    num_ = 0;
    gmx_bool bOK = 0;
    int ok = read_first_xtc(xtc, &numAtoms_, &step_, &time_, box_, &x_, &prec_, &bOK);
    // recentre on first atom
    recentreBox(0);

    std::ifstream gro(groname);
    if(!gro.is_open()){
        cout << "Error opening GRO file" << endl;
        exit(-1);
    }

    // first line of GRO is the system name
    char line[40];
    gro.getline(line, 40);
    cout << line << endl;
    // second line is the number of atoms
    int gro_num_atoms;
    gro >> gro_num_atoms;

    if(gro_num_atoms != numAtoms_){
        cout << "XTC num atoms:" << numAtoms_ << endl;
        cout << "GRO num atoms:" << gro_num_atoms << endl;
        cout << "Number of atoms declared in XTC file "
                "is not the same as declared in GRO file" << endl;
        throw std::runtime_error("Number of atoms does not match");
    }else{
        cout << "Found " << numAtoms_ << " atoms" << endl;
    }

    // print box vectors
    cout << "Box vectors" << endl;
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            printf("%8.4f", box_[i][j]);
        }
        cout << endl;
    }

    int res_interesting = 0;
    Parser parser(cfgname);
    vector<string> parse_buffer;
    parser.getLineFromSection("residues", &parse_buffer);
    res_interesting = std::stoi(parse_buffer[0]);
    cout << "Mapping first " << res_interesting << " residues" << endl;

    string res_name_new="", res_name_last="";
    int res_loc = -1;
    int res_num_atoms = 0;
    atoms_.resize(numAtoms_);
    // setup done, start reading in atoms from GRO
    for(int i = 0; i < numAtoms_; i++){
        atoms_[i] = Atom(i);
        string tmp_atom_type;
        try{
            gro >> res_name_new >> tmp_atom_type >> atoms_[i].atom_num;
        }catch(std::exception &e){
            cout << "Read from GRO file failed" << endl;
            exit(-1);
        }
        int tmp_int;
        sscanf(res_name_new.c_str(), "%d", &tmp_int);
        atoms_[i].atom_type = std::to_string(tmp_int) + tmp_atom_type;
        // skip rest of line
        gro.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        atoms_[i].atom_num--;
        atoms_[i].resid = res_name_new;

        // found a new residue
        if(res_name_new != res_name_last){
            res_loc++;
            residues_.push_back(Residue(res_name_new));
            res_num_atoms = 0;

            //TODO what if the residues we want aren't at the beginning
            if(res_loc < res_interesting + 1 && res_loc != 0){
                //TODO tidy up these - I want to print them, but nicely
                residues_[res_loc - 1].num_atoms = res_num_atoms;
                numAtomsTrack_ += residues_[res_loc - 1].atoms.size();
            }
        }

        residues_[res_loc].atoms.push_back(atoms_[i].atom_num);
        residues_[res_loc].atom_names.push_back(atoms_[i].atom_type);
        nameToNum_.emplace(atoms_[i].atom_type, i);
        numToName_.emplace(i, atoms_[i].atom_type);
        memcpy(atoms_[i].coords, x_[i], 3 * sizeof(float));

        // default values
        atoms_[i].charge = 0.f;
        atoms_[i].mass = 1.f;

        res_name_last = res_name_new;
        res_num_atoms++;
    }
    cout << "Mapping first " << numAtomsTrack_ << " atoms" << endl;

    // Process topology file
    vector<string> substrs;
    Parser top_parser(topname);
    for(int i=0; i<numAtomsTrack_; i++){
        // read data from topology file for each atom we care about (not solvent)
        // check that we're reading the atoms in the same order
        // internal atom name is the res # and atom name from top/gro
        top_parser.getLineFromSection("atoms", &substrs);
        string tmp_string = substrs[2] + substrs[4];
        assert(tmp_string == atoms_[i].atom_type);
        atoms_[i].charge = float(atof(substrs[6].c_str()));
        atoms_[i].mass = float(atof(substrs[7].c_str()));
    }
    gro.close();

    if(ok && bOK) isSetup_ = true;
    printAtoms(numAtomsTrack_);
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
    for(int i = 0; i < numAtoms_; i++){
        // overwrite coords of atoms stored in the current Frame
        memcpy(atoms_[i].coords, x_[i], 3 * sizeof(float));
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

void Frame::printAtoms(const int n){
    assert(isSetup_);
    cout.setf(std::ios::fixed);
    cout.precision(4);
    int i = 0;
    cout << "Name\tMass\tChrg\tPosx\tPosy\tPosz" << endl;
    for(Atom &atom : atoms_){
        cout << atom.atom_type << "\t" << atom.mass << "\t" << atom.charge << "\t";
        cout << atom.coords[0] << "\t" << atom.coords[1] << "\t" << atom.coords[2];
        cout << endl;
        i++;
        if(n > 0 && i >= n) break;
    }
}
