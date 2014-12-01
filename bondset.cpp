#include "bondset.h"

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

using std::vector;
using std::string;
using std::ifstream;

vector<string> &split(const string &s, char delim, vector<string> &elems);

vector<string> split(const string &s, char delim);

BondSet::BondSet(){
}

bool BondSet::fromFile(string filename){
    throw std::logic_error("Not implemented");
    bool ok = false;
    ifstream file;
    file.open(filename);
    if(file.is_open()){

    }
    return ok;
}

//TODO fix this; example from SE that doesn't work?
vector<string> &split(const string &s, char delim, vector<string> &elems){
    std::stringstream ss(s);
    string item;
    while(std::getline(ss, item, delim)){
        elems.push_back(item);
    }
    return elems;
}


vector<string> split(const std::string &s, char delim){
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}