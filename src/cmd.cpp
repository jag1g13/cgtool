#include "cmd.h"

#include <iostream>
#include <vector>

#include <boost/algorithm/string.hpp>

using std::string;
using std::cout;
using std::endl;
using std::vector;
using boost::algorithm::trim;

namespace po = boost::program_options;

CMD::CMD(const string &help_header, const string &help_string, const int argc, const char *argv[]){
    helpString_ = help_string;
    vector<string> lines;
    vector<string> parts(3);

    desc_.add_options()("help", "show this help text");
    boost::split(lines, helpString_, boost::is_any_of("\n"));
    for(const string &line : lines){
        boost::split(parts, line, boost::is_any_of("\t"));
        const string arg = boost::trim_left_copy_if(parts[0], boost::is_any_of("-"));
        desc_.add_options()(arg.c_str(), po::value<string>(), parts[1].c_str());
        stringArgs_.emplace(arg, parts[2]);
    }

    po::store(po::parse_command_line(argc, argv, desc_), options_);
    po::notify(options_);
    if (options_.count("help")) {
        cout << help_header;
        cout << desc_ << "\n";
        exit(0);
    }
}

const std::string CMD::getStringArg(const std::string &arg){
    if(options_.count(arg)){
        if(options_.count("dir")){
            return options_["dir"].as<string>() + "/" + options_[arg].as<string>();
        }else{
            return options_[arg].as<string>();
        }
    }
    if(stringArgs_.count(arg)){
        if(options_.count("dir")){
            return options_["dir"].as<string>() + "/" + stringArgs_[arg];
        }else{
            return stringArgs_[arg];
        }
    }
    cout << "Necessary argument not provied" << endl;
    exit(0);
}

const bool CMD::getBoolArg(const std::string &arg){
    // If it was set true by the user return it, otherwise return default value
    if(options_.count(arg)) return true;
    return bool(boolArgs_.count(arg));
}
