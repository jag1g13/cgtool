#include "cg_map.h"

#include <fstream>
#include <iostream>     // only needed for testing, get rid of this when done

//#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

#define DEBUG true

using std::cout;
using std::endl;

//vector<vector<string>> tokenize_file(string filename);

CGMap::CGMap(){
}

CGMap::CGMap(string filename){
    fromFile(filename);
}

bool CGMap::fromFile(string filename){
    bool status = 1;
    std::ifstream map_file(filename);
    string line;
    vector<string> substrs;
    if(!map_file.is_open()) return 0;   // couldn't open file
    while(getline(map_file, line)){
        if(line[0] == ';' || line[0] == '#') continue;  // skip comments
        if(line == "") continue;                        // line is empty, ignore it
        boost::split(substrs, line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
        BeadMap new_bead;
        new_bead.cg_bead = substrs[0];
        new_bead.atoms = vector<string>(substrs.begin() + 1, substrs.end());
        mapping.push_back(new_bead);
    }
    num_beads = mapping.size();
    if(DEBUG){
        for(auto &i : mapping){
            std::cout << i.cg_bead << " contains";
            for(auto &j : i.atoms){
                std::cout << " " << j;
            }
            std::cout << std::endl;
        }
    }
    return status;
}

void CGMap::initFrame(const Frame *aa_frame, Frame *cg_frame){
    cg_frame->allocateAtoms(num_beads);
    for(std::vector<BeadMap>::iterator i = mapping.begin(); i != mapping.end(); ++i){
        for(std::vector<string>::iterator j = i->atoms.begin(); j != i->atoms.end(); ++j){
            if(DEBUG) cout << *j << endl;
            //i->atom_nums.push_back();
        }
    }
}

bool CGMap::apply(const Frame *aa_frame, Frame *cg_frame){
    throw std::logic_error("Not implemented");

    bool status = true;
    return status;
}


//vector<vector<string>> tokenize_file(string filename){
//    vector<vector<string>> result;
//    string line;
//    std::ifstream map_file;
//    map_file.open(filename);
//    if(map_file.is_open()){
//        while(map_file >> line){
//            if(line[0] == ';' || line[0] == '#') continue;  // skip comments
//            boost::tokenizer<> tokens(line);                // otherwise split into tokens
//            for(boost::tokenizer<>::iterator token=tokens.begin(); token!=tokens.end(); ++token){
//
//            }
//        }
//    }
//}
