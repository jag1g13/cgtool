#ifndef PARSER_H_
#define PARSER_H_

#include <fstream>
#include <string>
#include <vector>

class Parser{
private:
public:
    /** The file currently being read */
    std::ifstream file_;
    /** The name of the file currently being read */
    std::string filename_;
    /** The line currently being read */
    std::string line_;
    /** The section of the file currently being read */
    std::string section_ = "";
    /** Is EOF? */
    bool eof_;

    /**
    * \brief Constructor for a blank Parser
    */
    Parser();

    /**
    * \brief Constructor for a Parser which will open a file and prepare for reading
    *
    * Raises an exception if file cannot be opened
    */
    Parser(std::string filename);

    /**
    * \brief Opens a file if the Parser was constructed without one
    *
    * Returns true if file has been opened successfully
    */
    bool openFile(std::string filename);

    /**
    * \brief Reads a line from file and splits it into tokens
    *
    * string section and Vector<string> will be filled with
    * the section of the file currently being read and the
    * tokens found on the current line
    */
    bool getLine(std::string *section, std::vector<std::string> *tokens);
};

#endif