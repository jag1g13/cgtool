#include "itp_writer.h"

#include <iostream>

using std::string;
using std::fprintf;
using std::cout;
using std::endl;
using std::vector;

ITPWriter::ITPWriter(string name){
    name_ = name;
    itp_ = std::fopen(name.c_str(), "w");
    if(itp_ == NULL){
        cout << "Could not open itp file for writing" << endl;
        exit(-1);
    }
    fprintf(itp_, "%s", header_.c_str());
}

void ITPWriter::newSection(std::string section_name){
    section_ = section_name;
    fprintf(itp_, "\n[ %s ]\n", section_.c_str());
}

void ITPWriter::printAtoms(const vector<BeadMap> &mapping){
    newSection("atoms");
    for(BeadMap bead : mapping){
        fprintf(itp_, "%6i %10s %6i %6s %6s %6i %10.4f %10.4f;\n",
                bead.num, bead.type.c_str(), bead.num, bead.name.c_str(),
                bead.name.c_str(), bead.num, bead.charge, bead.mass);
    }
}

//void ITPWriter::printBonds(const vector<> &bonds)