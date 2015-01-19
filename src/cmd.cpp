#include "cmd.h"

#include <iostream>
#include <sstream>

#include <boost/algorithm/string.hpp>

using std::string;
using std::cout;
using std::endl;
using boost::algorithm::trim;

/*
* Naval Fate.
*
* Usage:
*   naval_fate ship new <name>...
*   naval_fate ship <name> move <x> <y> [--speed=<kn>]
*   naval_fate ship shoot <x> <y>
*   naval_fate mine (set|remove) <x> <y> [--moored|--drifting]
*   naval_fate -h | --help
*   naval_fate --version
*
* Options:
*   -h --help     Show this screen.
*   --version     Show version.
*   --speed=<kn>  Speed in knots [default: 10].
*   --moored      Moored (anchored) mine.
*   --drifting    Drifting mine.
 */

CMD::CMD(string help_string){
    std::istringstream stream(help_string);
    while(!stream.eof()){
        string line;
        getline(stream, line);
        trim(line);
    }
}

CMD::CMD(){

}

bool CMD::parseArguments(const int argc, const char *argv[]){
//    throw std::runtime_error("Not implemented");
    string current_arg;
    for(int i=1; i<argc; i++){
        if(argv[i][0] == '-'){
            // new argument
            if(argv[i][1] == '-'){
                current_arg = argv[i][2];
            }else{
            current_arg = argv[i][1];
            }
            if(current_arg == "h" || current_arg == "help") help();
            if(argTypes_[current_arg] == ArgType::BOOL) boolArgs_[current_arg] == true;
            cout << current_arg << endl;

        }else{
            // if not new arg - must be value
            switch(argTypes_[current_arg]){
                case ArgType::STRING:
                    stringArgs_[current_arg] = argv[i];
                    break;
                case ArgType::INT:
                    intArgs_[current_arg] = std::stoi(argv[i]);
                    break;
                case ArgType::FLOAT:
                    floatArgs_[current_arg] = std::stof(argv[i]);
                    break;
                case ArgType::BOOL:
                    // do nothing, already dealt with
                    break;
            }
        }
    }
    return true;
}

void CMD::help(){
    cout << helpString_ << endl;
    exit(0);
}

const std::string CMD::getStringArg(const std::string arg){
    throw std::runtime_error("Not implemented");

}

const int CMD::getIntArg(const std::string arg){
    throw std::runtime_error("Not implemented");

}

const float CMD::getFloatArg(const std::string arg){
    throw std::runtime_error("Not implemented");

}

const bool CMD::getBoolArg(const std::string arg){
    throw std::runtime_error("Not implemented");

}
