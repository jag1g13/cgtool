#ifndef CMD_H_
#define CMD_H_

#include <string>

/**
* \brief Object to handle input to programs from the command line
* Parses the help text of the program to extract command line arguments.
* Parses the command line input of the program and stores all arguments into a dictionary.
* The program can query these dictionaries for the presence and values of arguments.
* Passing unknown arguments to the program will NOT cause it to fail.
*/
static class CMD{
protected:
    std::map<std::string, std::string> stringArgs_;
    std::map<std::string, int> intArgs_;
    std::map<std::string, float> floatArgs_;
    std::map<std::string, bool> boolArgs_;
public:
    /** \brief Constructor to parse the program help text */
    CMD(std::string help_string);

    /** \brief Parses arguments from the command line input */
    bool parseArguments(const int argc, const char* argv[]);

    /** \brief Return the value of named argument.
    * If argument was not provided by the user the default value will be used.
    * If there is not default value, print an error
    */
    const std::string getStringArg(const std::string arg);

    /** \brief Return the value of named argument.
    * If argument was not provided by the user the default value will be used.
    * If there is not default value, print an error
    */
    const int getIntArg(const std::string arg);

    /** \brief Return the value of named argument.
    * If argument was not provided by the user the default value will be used.
    * If there is not default value, print an error
    */
    const float getFloatArg(const std::string arg);

    /** \brief Return the value of named argument.
    * If argument was not provided by the user the default value will be used.
    * If there is not default value, print an error
    */
    const bool getBoolArg(const std::string arg);
};

#endif