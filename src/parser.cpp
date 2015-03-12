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

bool Parser::getLine(string &section, vector <string> &tokens){
    while(true){
        eof_ = !getline(file_, line_);
        boost::trim(line_);
        if(eof_) return !eof_;                              // stop if we hit eof
        if(line_[0] == ';' || line_[0] == '#') continue;    // skip comments
        if(line_ == "") continue;                           // line is empty, ignore it
        switch(format_){
            case ParserFormat::GROMACS:
                if(line_[0] == '['){                        // line is a section header
                    section_ = line_.substr(line_.find_first_of('[')+1, line_.find_last_of(']')-1);
                    boost::trim(section_);
                    continue;
                }
                break;

            case ParserFormat::LAMMPS:
                throw std::runtime_error("Not implemented");
            };
        break;                                              // line isn't empty, accept it
    }
    boost::split(tokens, line_, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
    for(string tok : tokens) boost::trim(tok);
    section = section_;
    return !eof_;       // return true if there is still file to read
}

bool Parser::findSection(const string find){
    string section = "";
    vector<string> token_buffer;
    while(section != find){
        if(!getLine(section, token_buffer)) return false;
    }
    return true;
}

bool Parser::getLineFromSection(const string find, vector<string> &tokens){
    //TODO don't read in anything if there isn't a line - why does it do this?
    string section_buffer;
    while(getLine(section_buffer, tokens)){
        if(section_ == find) return true;
    }
    rewind();
    return false;
}

void Parser::rewind(){
    // clear eof and rewind
    eof_ = false;
    file_.clear();
    file_.seekg(0, std::ios::beg);
}
