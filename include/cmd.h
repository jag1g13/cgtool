#ifndef CMD_H_
#define CMD_H_

#include "cmd_abstract.h"

#include <boost/program_options.hpp>

/**
* \brief Object to handle input to programs from the command line.
*/
class CMD : public CMD_Abstract{
protected:
    /** Store options from Boost program_options */
    boost::program_options::variables_map options_;
    boost::program_options::options_description desc_;

public:
    /** \brief Constructor to parse the program help text */
    CMD(const std::string &help_header, const std::string &help_string,
        const std::string &compile_info, const int argc, const char *argv[]);

    /** \brief Empty constructor, does nothing */
    CMD(){};

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
