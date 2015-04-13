#include "cmd.h"

#include <iostream>
#include <vector>

#include <sysexits.h>

#include <boost/algorithm/string.hpp>

using std::string;
using std::cout;
using std::endl;
using std::vector;
using boost::algorithm::trim;
using std::stoi;

namespace po = boost::program_options;

CMD::CMD(const string &help_header, const string &help_string, const int argc, const char *argv[]){
    helpString_ = help_string;
    vector<string> lines;
    vector<string> parts(3);

    // Split help string and parse it into options and default values
    desc_.add_options()("help,h", "show this help text");
    boost::split(lines, helpString_, boost::is_any_of("\n"));
    for(const string &line : lines){
        boost::split(parts, line, boost::is_any_of("\t"));
        const string arg = boost::trim_left_copy_if(parts[0], boost::is_any_of("-"));
        switch(static_cast<ArgType>(stoi(parts[3]))){
            case ArgType::PATH:
            case ArgType::STRING:
                desc_.add_options()((arg+","+arg[0]).c_str(),
                                    po::value<string>()->default_value(parts[2].c_str()),
                                    parts[1].c_str());
                stringArgs_.emplace(arg, parts[2]);
                break;
            case ArgType::INT:
                desc_.add_options()((arg+","+arg[0]).c_str(),
                                    po::value<int>()->default_value(stoi(parts[2])),
                                    parts[1].c_str());
                intArgs_.emplace(arg, stoi(parts[2]));
                break;
            case ArgType::FLOAT:
                break;
            case ArgType::BOOL:
                break;
        }
    }

    try{
        po::store(po::parse_command_line(argc, argv, desc_), options_);
    }catch(po::error e){
        cout << "Unrecognised command line argument\n" << endl;
        exit(EX_USAGE);
    }

    //command_line_parser(argc, argv).options(desc).allow_unregistered().run();
    po::notify(options_);
    if (options_.count("help")) {
        cout << help_header << endl;
        cout << "Arguments:" << endl;
        cout << desc_ << endl;
        exit(EX_OK);
    }
    if(options_.count("dir")) cout << "Files in: " << options_["dir"].as<string>() << endl;
}

const std::string CMD::getFileArg(const string &arg){
    // Was the argument passed in from the command line?
    if(options_.count(arg)){
        // All strings are file paths - append <dir> if user gave it
        if(options_.count("dir")){
            return options_["dir"].as<string>() + "/" + options_[arg].as<string>();
        }else{
            return options_[arg].as<string>();
        }
    }

    // Does the argument have a default value?
    if(stringArgs_.count(arg)){
        if(options_.count("dir")){
            return options_["dir"].as<string>() + "/" + stringArgs_[arg];
        }else{
            return stringArgs_[arg];
        }
    }

    // Can't find it - error
    cout << "Necessary argument not provied" << endl;
    exit(EX_USAGE);
}

const bool CMD::getBoolArg(const string &arg){
    // If it was set true by the user return it, otherwise return default value
    if(options_.count(arg)) return options_[arg].as<bool>();
    return boolArgs_[arg];
}

const int CMD::getIntArg(const string &arg){
    // If it was set by the user return it, otherwise return default value
    cout << options_.count(arg) << endl;
    cout << options_[arg].as<int>() << endl;
    if(options_.count(arg)) return options_[arg].as<int>();
    cout << "Returning from getIntArg" << cout;
    return 0;
}
