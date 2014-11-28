#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include "bondset.h"

using std::vector;
using std::string;
using std::ifstream;

vector<string> &split(const string &s, char delim, vector<string> &elems);
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
vector<string> &split(const string &s, char delim, vector<string> &elems) {
    std::stringstream ss(s);
    string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


vector<string> split(const std::string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}