//
// Created by james on 28/05/15.
//

#ifndef CGTOOL_COMMON_H
#define CGTOOL_COMMON_H

#include <string>
#include <vector>
#include <map>

#include "residue.h"
#include "frame.h"
#include "cg_map.h"
#include "itp_writer.h"
#include "parser.h"
#include "field_map.h"
#include "cmd.h"
#include "rdf.h"

struct CheckedFile{
    std::string name = "";
    bool exists = false;
};

struct DoFunction{
    bool on = false;
    int freq = 1;
    std::map<std::string, int> intProperty;
    std::map<std::string, double> doubleProperty;
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
    double lastUpdate_;

    // Input files
    std::map<std::string, CheckedFile> inputFiles_;

    // Run control
    int currFrame_ = 1;
    int numFramesMax_ = -1;
    int wholeXTCFrames_ = -1;
    bool untilEnd_ = true;
    std::map<std::string, DoFunction> doFunction_;

    // Objects
    std::vector<Residue> residues_;
    Frame    *frame_;
    Frame    *cgFrame_;
    BondSet  *bondSet_;
    CGMap    *cgMap_;
    RDF      *rdf_;
    FieldMap *field_;

    // Progress updates
    const int updateFreq_[10] = {1, 2, 5, 10, 20, 50, 100, 200, 500, 1000};
    int updateLoc_ = 0;

    // Function to perform in main loop
//    std::function<void(void)> *mainLoop;



public:
    Common();

    void setHelpStrings(const std::string &version, const std::string &header,
                        const std::string &options);

    void collectInput(const int argc, const char *argv[],
                      const std::vector<std::string> &req_files,
                      const std::vector<std::string> &opt_files);

    void findDoFunctions();

    void getResidues();

    void setupObjects();

    void doMainLoop();

    void updateProgress();

    void mainLoop();

    int do_stuff(const int argc, const char *argv[]);
};

#endif //CGTOOL_COMMON_H
