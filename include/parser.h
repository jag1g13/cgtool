#ifndef PARSER_H_
#define PARSER_H_

#include <fstream>
#include <string>
#include <vector>

enum class ParserFormat{GROMACS, LAMMPS};

/**
* \brief Parses input files for comments, section headers and data lines
*/
class Parser{
private:
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
    /** Expected file format - GROMACS style or LAMMPS style */
    ParserFormat format_ = ParserFormat::GROMACS;

    /**
    * \brief Reads a line from file and splits it into tokens
    *
    * Will skip over empty lines and comments and read section headers transparently
    * Parses the next data line and fills a vector<string> of tokens
    */
    bool getLine(std::string &section, std::vector<std::string> &tokens);

    /** \brief Rewind to start of file */
    void rewind();

public:
    /**
    * \brief Constructor for a Parser which will open a file and prepare for reading
    * \throws runtime_error if file cannot be opened
    */
    Parser(const std::string filename, const ParserFormat format=ParserFormat::GROMACS);

    /** \brief Search through a config file for a particular section
    * Returns false if section cannot be found
    */
    bool findSection(const std::string find);

    /**\brief Search through config file for a particular section and pass back lines
    * Once it reaches the end of the file, it rewinds to the beginning and returns */
    bool getLineFromSection(const std::string find, std::vector<std::string> &tokens);
};

#endif