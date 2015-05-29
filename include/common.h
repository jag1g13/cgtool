//
// Created by james on 28/05/15.
//

#ifndef CGTOOL_COMMON_H
#define CGTOOL_COMMON_H

#include <string>
#include <vector>
#include <map>

struct CheckedFile{
    std::string name = "";
    bool exists = false;
};

struct DoFunction{
    bool on = false;
    int freq = 1;
};

class Common{
protected:
    // Help texts
    std::string versionString_;
    std::string helpHeader_;
    std::string helpOptions_;

    // Timing
    double veryStart_;
    double sectionStart_;

    // Input files
    std::map<std::string, CheckedFile> inputFiles_;

    // Run control
    int numFramesMax_ = -1;
    std::map<std::string, DoFunction> doFunction_;

    // Objects
    std::vector<Residue> residues_;

public:
    Common();

    void setHelpStrings(const std::string &version, const std::string &header,
                        const std::string &options);

    void collectInput(const int argc, const char *argv[],
                      const std::vector<std::string> &req_files,
                      const std::vector<std::string> &opt_files);

    void findDoFunctions();

    void getResidues();

    int do_stuff(const int argc, const char *argv[]);
};

#endif //CGTOOL_COMMON_H
