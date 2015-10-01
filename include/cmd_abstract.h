#ifndef CMD_ABSTRACT_H_
#define CMD_ABSTRACT_H_

#include <string>

enum class ArgType{PATH, STRING, INT, FLOAT, BOOL};

/**
* \brief Object to handle input to programs from the command line.
*/
class CMD_Abstract{
protected:
    /** Program help string.  Should be parsed to generate arguments */
    std::string helpString_;

public:
    /** \brief Return the value of named filepath argument.
    * If argument was not provided by the user the default value will be used.
    * If there is no default value, print an error
    */
    virtual std::string getStringArg(const std::string &arg) const = 0;

    /** \brief Return the value of a named boolean argument.
    * If argument was not provided by the user the default value will be used.
    * If there is no default value, print an error
    */
    virtual bool getBoolArg(const std::string &arg) const = 0;

    /** \brief Return the value of a named integer argument.
    * If argument was not provided by the user the default value will be used.
    * If there is no default value, print an error
    */
    virtual int getIntArg(const std::string &arg) const = 0;
};

#endif
