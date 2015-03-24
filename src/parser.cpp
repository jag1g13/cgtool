#include "parser.h"

#include <iostream>

#include <boost/algorithm/string.hpp>

using std::string;
using std::vector;
using std::cout;
using std::endl;

Parser::Parser(const string filename, const ParserFormat format) {
    //TODO preprocess file to include ITPs
    format_ = format;
    filename_ = filename;
    file_.open(filename);
    if (!file_.is_open()) throw std::runtime_error("File " + filename + " could not be opened");
}

Parser::~Parser(){
    if(file_.is_open()) file_.close();
}

bool Parser::getLine(vector <string> &tokens){
    while(true){
        eof_ = !getline(file_, line_);

        // Stop if we hit eof
        if(eof_) return false;
        boost::trim(line_);

        // Skip comments
        if(line_[0] == ';' || line_[0] == '#') continue;

        // Line is empty, ignore it
        if(line_ == "") continue;

        switch(format_){
            case ParserFormat::GROMACS:
                // Line is a section header
                if(line_[0] == '['){
                    section_ = line_.substr(line_.find_first_of('[')+1, line_.find_last_of(']')-1);
                    boost::trim(section_);
                    continue;
                }
                break;

            case ParserFormat::LAMMPS:
                throw std::logic_error("Not implemented");
            };

        // Line isn't empty, accept it
        break;
    }

    // Separate and trim whitespace from tokens
    boost::split(tokens, line_, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
    for(string &tok : tokens) boost::trim(tok);

    // Return true if there is still file to read
    return true;
}

bool Parser::findSection(const string find){
    vector<string> token_buffer;
    while(section_ != find){
        if(!getLine(token_buffer)) return false;
    }
    return true;
}

bool Parser::getLineFromSection(const string find, vector<string> &tokens){
    // We're looking for a new section - it might be above the last one
    if(find != findPrevious_) rewind();
    findPrevious_ = find;
    while(getLine(tokens)){
        if(section_ == find) return true;
    }
    rewind();
    return false;
}

void Parser::rewind(){
    // Clear eof and rewind
    eof_ = false;
    file_.clear();
    file_.seekg(0, std::ios::beg);
}
