#include "parser.h"

#include <iostream>

#include <boost/algorithm/string.hpp>

using std::string;
using std::vector;
using std::cout;
using std::endl;

Parser::Parser() {

}

Parser::Parser(string filename, ParserFormat format) {
    format_ = format;
    filename_ = filename;
    if (!openFile(filename)) throw std::runtime_error("File " + filename + " could not be opened");
}

bool Parser::openFile(string filename){
    //TODO preprocess file to include ITPs
    filename_ = filename;
    file_.open(filename);
    return file_.is_open();
}

bool Parser::getLine(string *section, vector <string> *tokens){
    while(true){
        eof_ = !getline(file_, line_);
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
    boost::trim(line_);
    boost::split(*tokens, line_, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
    for(string tok : *tokens) boost::trim(tok);
    *section = section_;
    return !eof_;       // return true if there is still file to read
}

bool Parser::findSection(const string find){
    string section = "";
    vector<string> token_buffer;
    while(section != find){
        if(!getLine(&section, &token_buffer)) return false;
    }
    return true;
}

bool Parser::findNextSection(){
    while(true){
        eof_ = !getline(file_, line_);
        if(eof_) return false;                              // stop if we hit eof
        switch(format_){
            case ParserFormat::GROMACS:
                if(line_[0] == '['){                        // line is a section header
                    section_ = line_.substr(line_.find_first_of('[')+1, line_.find_last_of(']')-1);
                    boost::trim(section_);
                    return true;
                }
                break;

            case ParserFormat::LAMMPS:
                throw std::runtime_error("Not implemented");
        }
    }
}

bool Parser::getLastLineFromSection(const string find, vector<string> *tokens){
    string section_buffer;
    while(getLine(&section_buffer, tokens)){
        if(section_ == find) break;
    }
    rewind();
    return false;
}

bool Parser::getLineFromSection(const string find, vector<string> *tokens){
    string section_buffer;
    while(getLine(&section_buffer, tokens)){
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
