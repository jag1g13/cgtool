#ifndef CMD_H_
#define CMD_H_

#include <string>
#include <map>
#include <boost/program_options.hpp>

enum class ArgType{PATH, STRING, INT, FLOAT, BOOL};

/**
* \brief Object to handle input to programs from the command line.
*
* Parses the help text of the program to extract command line arguments.
* Parses the command line input of the program and stores all arguments into a dictionary.
* The program can query these dictionaries for the presence and values of arguments.
* Passing unknown arguments to the program will NOT cause it to fail.
*
*/
class CMD{
protected:
    /** \brief Program help string.  Should be parsed to generate arguments */
    std::string helpString_;

    /** \brief Store options from Boost program_options */
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
    const std::string getFileArg(const std::string &arg);

    /** \brief Return the value of a named boolean argument.
    * If argument was not provided by the user the default value will be used.
    * If there is no default value, print an error
    */
    const bool getBoolArg(const std::string &arg);

    /** \brief Return the value of a named integer argument.
    * If argument was not provided by the user the default value will be used.
    * If there is no default value, print an error
    */
    const int getIntArg(const std::string &arg);
};

#endif