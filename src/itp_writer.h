#ifndef ITP_WRITER_H_
#define ITP_WRITER_H_

#include <fstream>
#include <string>

class ITPWriter{
protected:
    std::ofstream itp_;
    std::string name_;
    std::string section_;

public:
    /** Create an ITP file and prepare to write */
    ITPWriter(std::string name);

    /** Create a new section in the ITP file */
    void newSection(std::string section_name);
};

#endif