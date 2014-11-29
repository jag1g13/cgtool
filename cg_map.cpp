#include <fstream>

#include <boost/tokenizer.hpp>

#include "cg_map.h"

vector<vector<string>> tokenize_file(string filename);

bool CGMap::from_file(string filename){
    bool status = 0;
    return status;
}

bool CGMap::apply(Frame* aa_frame, Frame* cg_frame){
    bool status = true;
    return status;
}

vector<vector<string>> tokenize_file(string filename){
    vector<vector<string>> result;
    string line;
    std::ifstream map_file;
    map_file.open(filename);
    if(map_file.is_open()){
        while(map_file >> line){
            if(line[0] == ';' || line[0] == '#') continue;  // skip comments
            boost::tokenizer<> tokens(line);                // otherwise split into tokens
            for(boost::tokenizer<>::iterator token=tokens.begin(); token!=tokens.end(); ++token){

            }
        }
    }
}