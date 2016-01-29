#include "cmd_simple.h"

#include <iostream>

#include <sysexits.h>

#include <boost/algorithm/string.hpp>

using std::string;
using std::cout;
using std::endl;
using std::vector;
using boost::algorithm::trim;
using std::stoi;

CMDSimple::CMDSimple(const string &help_header, const string &help_string,
         const string &compile_info, const int argc, const char *argv[]){
    helpString_ = help_string;
    vector<string> lines;
    vector<string> parts(3);

    // Split help string and parse it into options and default values
    boost::split(lines, helpString_, boost::is_any_of("\n"));
    for(const string &line : lines){
        boost::split(parts, line, boost::is_any_of("\t"));
        const string arg = boost::trim_left_copy_if(parts[0], boost::is_any_of("-"));

        type_[arg] = static_cast<ArgType>(stoi(parts[2]));

        if(type_[arg] == ArgType::STRING){
            // String gets a short form
            shortForm_[arg] = arg[0];
            options_[arg] = "";
        }else{
            // Everything else gets a default value
            options_[arg] = parts[3];
        }

    }

    // Parse arguments
    for(int i=1; i<argc; i++){
        const string arg = boost::trim_left_copy_if(argv[i], boost::is_any_of("-"));
        if(options_.count(arg)){
            if(type_[arg] == ArgType::BOOL){
                options_[arg] = "true";
            }else{
                options_[arg] = argv[i+1];
            }

        }else{
            cout << "Unrecognised command line argument\n" << endl;
            cout << "Arguments:" << endl;
            exit(EX_USAGE);
        }
    }

    if(options_.count("help") || argc < 2){
        cout << help_header << endl;
        cout << "Arguments:" << endl;
        exit(EX_OK);
    }

    if(options_.count("version")){
        cout << compile_info << endl;
        exit(EX_OK);
    }
}

std::string CMDSimple::getStringArg(const string &arg) const{
    if(options_.count(arg)) return options_.at(arg);
    // No value or default - return empty string
    return "";
}

bool CMDSimple::getBoolArg(const string &arg) const{
    if(options_.count(arg)) return (options_.at(arg) != "");
    // No value or default - return false
    return false;
}

int CMDSimple::getIntArg(const string &arg) const{
    if(options_.count(arg)) return std::stoi(options_.at(arg));
    // No value or default - return 0
    return 0;
}
