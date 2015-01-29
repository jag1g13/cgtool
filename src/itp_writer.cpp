#include "itp_writer.h"

using std::string;

ITPWriter::ITPWriter(string name){
    name_ = name;
}

void ITPWriter::newSection(std::string section_name){
    section_ = section_name;
}