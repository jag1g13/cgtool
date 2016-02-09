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

    // Add help and version manually
    type_["help"] = ArgType::BOOL;
    shortForm_['h'] = "help";
    options_["help"] = "";
    type_["version"] = ArgType::BOOL;
    shortForm_['v'] = "version";
    options_["version"] = "";

    // Split help string and parse it into options and default values
    boost::split(lines, helpString_, boost::is_any_of("\n"));
    for(const string &line : lines){
        boost::split(parts, line, boost::is_any_of("\t"));
        string arg = parts[0];
        while(arg[0] == '-') arg.replace(0, 1, "");

        type_[arg] = static_cast<ArgType>(stoi(parts[2]));

        if(type_[arg] == ArgType::STRING){
            // String gets a short form
            shortForm_[arg.at(0)] = arg;
            options_[arg] = "";
        }else{
            // Everything else gets a default value
            options_[arg] = parts[3];
        }

    }

    // Parse arguments
    for(int i=1; i<argc; i++){
        string arg = argv[i];
        // Strip leading hyphen and check if short form
        if(arg[0] == '-'){
            arg.replace(0, 1, "");
            if(arg[0] == '-'){
                arg.replace(0, 1, "");
            }else{
                if(arg.length() == 1 && shortForm_.count(arg.at(0))){
                    arg = shortForm_[arg.at(0)];
                }
            }
        }

        // Register arguments in map
        if(options_.count(arg)){
            if(type_[arg] == ArgType::BOOL){
                options_[arg] = "true";
            }else{
                if(argc < i+2){
                    cout << "Argument parameter not given\n" << endl;
                    cout << help_header << endl;
                    exit(EX_USAGE);
                }
                options_[arg] = argv[i+1];
                i++;
            }

        }else{
            cout << "Unrecognised command line argument\n" << endl;
            cout << help_header << endl;
            exit(EX_USAGE);
        }
    }

    if(options_.at("help") != "" || argc < 2){
        cout << help_header << endl;
        exit(EX_OK);
    }

    if(options_.at("version") != ""){
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
