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
//    desc_.add_options()("help,h", "Show this help text");
//    desc_.add_options()("version,v", "Show version information");

    boost::split(lines, helpString_, boost::is_any_of("\n"));
    for(const string &line : lines){
        boost::split(parts, line, boost::is_any_of("\t"));
        const string arg = boost::trim_left_copy_if(parts[0], boost::is_any_of("-"));

        options_[arg] = NULL;
        type_[arg] = static_cast<ArgType>(stoi(parts[2]));

        if(type_[arg] == ArgType::STRING){
            // String gets a short form
            shortForm_[arg] = arg[0];
        }else{
            // Everything else gets a default value
            //TODO Probably need to sort out type here
            options_[arg] = parts(3);
        }

    }

    // Parse arguments
    for(int i=1; i<argc; i++){
        const string arg = boost::trim_left_copy_if(argv[i], boost::is_any_of("-"));
        if(options_.count(arg)){
            if(type_[arg] == ArgType::BOOL){
                options_[arg] = true;
            }else{
                optios_[arg] = argv[i+1];
            }

        }else{
            cout << "Unrecognised command line argument\n" << endl;
            cout << "Arguments:" << endl;
//        cout << desc_ << endl;
            exit(EX_USAGE);
        }
    }

//    po::notify(options_);
    if(options_.count("help") || argc < 2){
        cout << help_header << endl;
        cout << "Arguments:" << endl;
//        cout << desc_ << endl;
        exit(EX_OK);
    }

    if(options_.count("version")){
        cout << compile_info << endl;
        exit(EX_OK);
    }
}

std::string CMDSimple::getStringArg(const string &arg) const{
    if(options_.count(arg)) return options_[arg].as<string>();
    // No value or default - return empty string
    return "";
}

bool CMDSimple::getBoolArg(const string &arg) const{
    if(options_.count(arg)) return options_[arg].as<bool>();
    // No value or default - return false
    return false;
}

int CMDSimple::getIntArg(const string &arg) const{
    if(options_.count(arg)) return options_[arg].as<int>();
    // No value or default - return -1
    return -1;
}
