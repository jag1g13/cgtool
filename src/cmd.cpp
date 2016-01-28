#include "cmd.h"

#include <iostream>

#include <sysexits.h>

#include <boost/algorithm/string.hpp>

using std::string;
using std::cout;
using std::endl;
using std::vector;
using boost::algorithm::trim;
using std::stoi;

namespace po = boost::program_options;

CMD::CMD(const string &help_header, const string &help_string,
         const string &compile_info, const int argc, const char *argv[]){
    helpString_ = help_string;
    vector<string> lines;
    vector<string> parts(3);

    // Split help string and parse it into options and default values
    desc_.add_options()("help,h", "Show this help text");
    desc_.add_options()("version,v", "Show version information");

    boost::split(lines, helpString_, boost::is_any_of("\n"));
    for(const string &line : lines){
        boost::split(parts, line, boost::is_any_of("\t"));
        const string arg = boost::trim_left_copy_if(parts[0], boost::is_any_of("-"));
        switch(static_cast<ArgType>(stoi(parts[2]))){
            case ArgType::STRING:
                // String arguments get a short form
                desc_.add_options()((arg + "," + arg[0]).c_str(),
                                    po::value<string>(),
                                    parts[1].c_str());
                break;
            case ArgType::INT:
                desc_.add_options()((arg).c_str(),
                                    po::value<int>()->default_value(stoi(parts[3])),
                                    parts[1].c_str());
                break;
            case ArgType::FLOAT:
                break;
            case ArgType::BOOL:
                desc_.add_options()((arg).c_str(),
                                    po::value<bool>()->default_value(stoi(parts[3])),
                                    parts[1].c_str());
                break;
        }
    }

    try{
        po::store(po::parse_command_line(argc, argv, desc_), options_);
    }catch(po::error e){
        cout << "Unrecognised command line argument\n" << endl;
        cout << "Arguments:" << endl;
        cout << desc_ << endl;
        exit(EX_USAGE);
    }

    po::notify(options_);
    if(options_.count("help") || argc == 1){
        cout << help_header << endl;
        cout << "Arguments:" << endl;
        cout << desc_ << endl;
        exit(EX_OK);
    }

    if(options_.count("version")){
        cout << compile_info << endl;
        exit(EX_OK);
    }
}

std::string CMD::getStringArg(const string &arg) const{
    if(options_.count(arg)) return options_[arg].as<string>();
    // No value or default - return empty string
    return "";
}

bool CMD::getBoolArg(const string &arg) const{
    if(options_.count(arg)) return options_[arg].as<bool>();
    // No value or default - return false
    return false;
}

int CMD::getIntArg(const string &arg) const{
    if(options_.count(arg)) return options_[arg].as<int>();
    // No value or default - return 0
    return 0;
}
