#include "frame.h"

#include <iostream>
#include <sstream>

#include <math.h>
#include <assert.h>
#include <sysexits.h>

#include "xdrfile_xtc.h"

#include "parser.h"
#include "small_functions.h"

using std::string;
using std::vector;
using std::endl;
using std::stoi;
using std::stof;
using std::printf;
using std::map;

Frame::Frame(const Frame &frame, vector<Residue> *residues){
    num_ = frame.num_;
    name_ = frame.name_;
    prec_ = frame.prec_;
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

    if(!initFromXTC(xtcname)){
        printf("ERROR: Something went wrong with reading XTC file\n");
        exit(EX_UNAVAILABLE);
    };
    createAtoms(numAtoms_);

    // Populate atoms_
    if(!initFromGRO(groname)){
        printf("ERROR: Something went wrong with reading GRO file\n");
        exit(EX_UNAVAILABLE);
    };

}

Frame::~Frame(){
    isSetup_ = false;
    if(xtcOutput_) xdrfile_close(xtcOutput_);
    if(xtcInput_) xdrfile_close(xtcInput_);
    if(x_) delete[] x_;
}

void Frame::setupOutput(string xtcname, string topname){
    char mode[2] = {'r', 'w'};
    if(xtcname == "") xtcname = (*residues_)[0].resname + ".xtc";
    if(topname == "") topname = (*residues_)[0].resname + ".top";
    backup_old_file(xtcname);
    backup_old_file(topname);
    if(!xtcOutput_) xtcOutput_ = xdrfile_open(xtcname.c_str(), &mode[1]);
    if(!xtcOutput_) throw std::runtime_error("Could not open XTC output");

    if(!x_) x_ = new rvec[numAtoms_];
    if(!x_) throw std::runtime_error("Couldn't allocate memory for XTC output");

    std::ofstream top(topname);
    if(!top.is_open()) throw std::runtime_error("Could not open output TOP file");
    top << "; Include forcefield parameters" << endl;
    top << "#include \"" << (*residues_)[0].resname << ".itp\"" << endl << endl;
    top << "[ system ]" << endl;
    top << (*residues_)[0].resname << endl << endl;
    top << "[ molecules ]" << endl;
    top << (*residues_)[0].resname << "\t\t" << (*residues_)[0].num_residues << endl;
    top.close();
    outputSetup_ = true;
}

bool Frame::writeToXtc(){
    if(!outputSetup_) throw std::logic_error("Output has not been setup");
    // need to put atomic coordinates back into x_
    // either it's a CG frame and x_ is empty, or it's atomistic but may have been recentred
    for(int i=0; i<(*residues_)[0].total_atoms; i++){
        x_[i][0] = float(atoms_[i].coords[0]);
        x_[i][1] = float(atoms_[i].coords[1]);
        x_[i][2] = float(atoms_[i].coords[2]);
    }
    return exdrOK == write_xtc(xtcOutput_, (*residues_)[0].total_atoms, step_, time_, box_, x_, prec_);
}

bool Frame::initFromXTC(const string &xtcname){
    // How many atoms?  Prepare Frame for reading
    int status = read_xtc_natoms(xtcname.c_str(), &numAtoms_);
    if(status != exdrOK) return false;

    xtcInput_ = xdrfile_open(xtcname.c_str(), "r");
    num_ = 0;

    // Init system from XTC file - libxdrfile library
    x_ = new rvec[numAtoms_];
    status = read_xtc(xtcInput_, numAtoms_, &step_, &time_, box_, x_, &prec_);

    // Check box vectors
    boxType_ = BoxType::CUBIC;
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            if(i!=j && box_[i][j] > 0){
                boxType_ = BoxType::TRICLINIC;
            }
        }
    }
    if(boxType_ != BoxType::CUBIC) printf("NOTE: Box is not cubic\n");

    if(status == exdrOK) isSetup_ = true;
    return isSetup_;
}

bool Frame::initFromGRO(const string &groname){
    // Require that atoms have been created
    assert(atomHas_.created);

    std::ifstream gro(groname);
    if(!gro.is_open()) return false;

    // Read top two lines of GRO - system name and number of atoms
    string tmp;
    int num_atoms = 0;
    getline(gro, name_);
    getline(gro, tmp);
    num_atoms = stoi(tmp);
    assert(num_atoms == numAtoms_);

    // The columns of a GRO file
    vector<GROLine> grolines(num_atoms);

    // Read in GRO file and get num of different resnames
    int num_residues = 1;
    getline(gro, tmp);
    grolines[0].populate(tmp);
    for(int i=1; i<numAtoms_; i++){
        getline(gro, tmp);
        grolines[i].populate(tmp);
        if(grolines[i].resname.compare(grolines[i-1].resname)) num_residues++;
    }
    gro.close();

    if(residues_->size() < num_residues){
        printf("ERROR: Found %'d residue(s) not listed in CFG\n",
               num_residues - static_cast<int>(residues_->size()));
        exit(EX_CONFIG);
    }

    int resnum = 0, atomnum = 0;
    num_residues = 1;
    Residue *res = &((*residues_)[0]);
    res->set_start(0);

    // Final pass to populate residues
    atoms_[0].resnum = grolines[0].resnum;
    atoms_[0].atom_name = grolines[0].atomname;
    for(int i=1; i<numAtoms_; i++){
        atoms_[i].resnum = grolines[i].resnum;
        atoms_[i].atom_name = grolines[i].atomname;
        atoms_[i].coords[0] = grolines[i].coords[0];
        atoms_[i].coords[1] = grolines[i].coords[1];
        atoms_[i].coords[2] = grolines[i].coords[2];

        if(grolines[i].resnum != grolines[i-1].resnum) resnum++;
        atomnum++;

        if(grolines[i].resname.compare(grolines[i-1].resname)){
            res->set_resname(grolines[i-1].resname);
            res->set_num_residues(resnum);
            res->set_total_atoms(atomnum);
            res->set_num_atoms(res->total_atoms / res->num_residues);
            res->calc_total();
            res->populated = true;

            num_residues++;
            res = &((*residues_)[num_residues-1]);
            res->set_start(i);
            resnum = 0;
            atomnum = 0;
        }
    }

    // Complete final residue
    resnum++;
    atomnum++;
    res->set_resname(grolines[numAtoms_-1].resname);
    res->set_num_residues(resnum);
    res->set_total_atoms(atomnum);
    res->set_num_atoms(res->total_atoms / res->num_residues);
    res->calc_total();
    res->populated = true;

    atomHas_.atom_name = true;
    atomHas_.resnum = true;
    atomHas_.coords = true;

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

void Frame::copyCoordsIntoAtoms(int natoms){
    if(natoms < 0) natoms = numAtoms_;
    assert(atoms_.size() >= natoms);
    for(int i=0; i<natoms; i++){
        atoms_[i].coords[0] = x_[i][0];
        atoms_[i].coords[1] = x_[i][1];
        atoms_[i].coords[2] = x_[i][2];
    }
    atomHas_.coords = true;
}

void Frame::createAtoms(int natoms){
    if(natoms < 0) natoms = numAtoms_;
    atoms_.resize(natoms);
    atomHas_.created = true;
}

void Frame::pbcAtom(int natoms){
    if(natoms < 0) natoms = numAtoms_;
    for(int i=0; i<natoms; i++){
        for(int j=0; j<3; j++){
            // For each coordinate wrap around into box
            while(x_[i][j] < 0.) x_[i][j] += box_[j][j];
            while(x_[i][j] > box_[j][j]) x_[i][j] -= box_[j][j];
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
    assert(xtcInput_ != nullptr);
    int status = read_xtc(xtcInput_, numAtoms_, &step_, &time_, box_, x_, &prec_);
    if(status != exdrOK) return false;
    num_++;

    pbcAtom();
    copyCoordsIntoAtoms();
    return true;
}

void Frame::printAtoms(int natoms) const{
    assert(isSetup_);
    if(natoms == -1) natoms = numAtoms_;
    printf("  Num Name    Mass  Charge    Posx    Posy    Posz\n");
    for(int i=0; i<natoms; i++){
        printf("%5i%5s%8.3f%8.3f%8.3f%8.3f%8.3f\n",
               i, atoms_[i].atom_name.c_str(),
               atoms_[i].mass, atoms_[i].charge,
               atoms_[i].coords[0], atoms_[i].coords[1], atoms_[i].coords[2]);
    }
}

void Frame::printGRO(string filename, int natoms) const{
    assert(isSetup_);
    if(filename == "") filename = (*residues_)[0].resname + ".gro";
    if(natoms == -1) natoms = numAtoms_;
    backup_old_file(filename);

    FILE *gro = std::fopen(filename.c_str(), "w");
    if(gro == nullptr){
        printf("ERROR: Could not open gro file for writing\n");
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
    stream << "Generated by CGTOOL : " << (*residues_)[0].resname << "\n" << natoms << "\n";
    fprintf(gro, "%s", stream.str().c_str());
    for(int i=0; i < natoms; i++){
        fprintf(gro, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
                1+(i/(*residues_)[0].num_atoms), (*residues_)[0].resname.c_str(),
                atoms_[i].atom_name.c_str(), i+1,
                atoms_[i].coords[0], atoms_[i].coords[1], atoms_[i].coords[2]);

        for(int j=0; j < 3; j++){
            min[j] = fmin(min[j], atoms_[i].coords[j]);
            max[j] = fmax(max[j], atoms_[i].coords[j]);
        }
    }
    fprintf(gro, "%10.5f%10.5f%10.5f\n", max[0]-min[0], max[1]-min[1], max[2]-min[2]);

    fclose(gro);
}

void Frame::printBox() const{
    // Print box vectors - assume cubic
    for(int i=0; i<3; i++){
        printf("%8.4f", box_[i][i]);
    }
    printf("\n");
}

