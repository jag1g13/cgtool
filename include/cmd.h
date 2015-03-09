#ifndef CMD_H_
#define CMD_H_

#include <string>
#include <map>
#include <boost/program_options.hpp>

enum class ArgType{STRING, INT, FLOAT, BOOL};

/**
* \brief Object to handle input to programs from the command line
* Parses the help text of the program to extract command line arguments.
* Parses the command line input of the program and stores all arguments into a dictionary.
* The program can query these dictionaries for the presence and values of arguments.
* Passing unknown arguments to the program will NOT cause it to fail.
*
*/
class CMD{
protected:
    /** Stores string argument names */
    std::map<std::string, std::string> stringArgs_;

    /** Stores boolean arguments accessed by argument name */
    std::map<std::string, bool> boolArgs_;

    /** Program help string.  Should be parsed to generate arguments */
    std::string helpString_;

    /** Store options from Boost program_options */
    boost::program_options::variables_map options_;
    boost::program_options::options_description desc_;

public:
    /** \brief Constructor to parse the program help text */
    CMD(const std::string &help_string, const int argc, const char *argv[]);

    /** \brief Empty constructor, does nothing */
    CMD(){};

    /** \brief Return the value of named argument.
    * If argument was not provided by the user the default value will be used.
    * If there is not default value, print an error
    */
    const std::string &getStringArg(const std::string &arg);

    /** \brief Return the value of named argument.
    * If argument was not provided by the user the default value will be used.
    * If there is not default value, print an error
    */
    const int getIntArg(const std::string &arg);

    /** \brief Return the value of named argument.
    * If argument was not provided by the user the default value will be used.
    * If there is not default value, print an error
    */
    const float getFloatArg(const std::string &arg);

    /** \brief Return the value of named argument.
    * If argument was not provided by the user the default value will be used.
    * If there is not default value, print an error
    */
    const bool getBoolArg(const std::string &arg);
};

#endif