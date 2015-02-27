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

CMD::CMD(const string &help_string, const int argc, const char *argv[]){
//    const string help_options =
//            "--xtc\tGROMACS xtc file\tmd.xtc\n"
//            "--itp\tGROMACS itp file\ttopol.top\n"
//            "--cfg\tCGTOOL mapping file\tcg.cfg";

    helpString_ = help_string;
    vector<string> lines;
    vector<string> parts(3);

//    desc_.options_description("Allowed options");
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
        cout << desc_ << "\n";
    }
    exit(0);
}

CMD::CMD(){

}

void CMD::help(){
    cout << helpString_ << endl;
    exit(0);
}

const std::string &CMD::getStringArg(const std::string &arg){
    throw std::runtime_error("Not implemented");

}

const int CMD::getIntArg(const std::string &arg){
    throw std::runtime_error("Not implemented");

}

const float CMD::getFloatArg(const std::string &arg){
    throw std::runtime_error("Not implemented");

}

const bool CMD::getBoolArg(const std::string &arg){
    throw std::runtime_error("Not implemented");

}
