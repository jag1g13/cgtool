#ifndef CMD_SIMPLE_H_
#define CMD_SIMPLE_H_

#include "cmd_abstract.h"

#include <map>

/**
* \brief Object to handle input to programs from the command line.
*/
class CMDSimple : public CMDAbstract{
protected:
    /** Store options */
    std::map<std::string, std::string> options_;
    std::map<std::string, ArgType> type_;
    std::map<std::string, std::string> shortForm_;
    std::map<std::string, std::string> descriptions_;

public:
    /** \brief Constructor to parse the program help text */
    CMDSimple(const std::string &help_header, const std::string &help_string,
        const std::string &compile_info, const int argc, const char *argv[]);

    /** \brief Empty destructor, does nothing */
    ~CMDSimple(){};

    /** \brief Return the value of named filepath argument.
    * If argument was not provided by the user the default value will be used.
    * If there is no default value, print an error
    */
    std::string getStringArg(const std::string &arg) const;

    /** \brief Return the value of a named boolean argument.
    * If argument was not provided by the user the default value will be used.
    * If there is no default value, print an error
    */
    bool getBoolArg(const std::string &arg) const;

    /** \brief Return the value of a named integer argument.
    * If argument was not provided by the user the default value will be used.
    * If there is no default value, print an error
    */
    int getIntArg(const std::string &arg) const;
};

#endif
