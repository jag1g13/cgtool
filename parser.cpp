#include "parser.h"

#include <boost/algorithm/string.hpp>

using std::string;
using std::vector;

Parser::Parser() {

}

Parser::Parser(string filename) {
    filename_ = filename;
    if (!openFile(filename)) throw std::runtime_error("File " + filename + " could not be opened");
}

bool Parser::openFile(string filename){
    file_.open(filename);
    return file_.is_open();
}

bool Parser::getLine(string *section, vector <string> *tokens){
    while(true){
        eof_ = !getline(file_, line_);
        if(eof_) return !eof_;                              // stop if we hit eof
        if(line_[0] == ';' || line_[0] == '#') continue;    // skip comments
        if(line_ == "") continue;                           // line is empty, ignore it
        if(line_[0] == '['){                                // line is a section header
            section_ = line_.substr(line_.find_first_of('['), line_.find_last_of(']'));
            continue;
        }
        break;                                              // line isn't empty, accept it
    }
    boost::split(*tokens, line_, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
    *section = section_;
    return !eof_;       // return true if there is still file to read
}