#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include "bondset.h"

using std::vector;
using std::string;
using std::ifstream;

vector<string> split(const string &s, char delim);

BondSet::BondSet(){
}

bool BondSet::from_file(string filename){
    ifstream file;
    file.open(filename);
    if(file.is_open()){

    }
}

//TODO fix this; example from SE that doesn't work?
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    //split(s, delim, elems);
    return elems;
}