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
#include "membrane.h"

struct CheckedFile{
    std::string name = "";
    bool exists = false;
};

struct DoFunction{
    bool on = false;
    int freq = 1;
    std::map<std::string, int> intProperty;
    std::map<std::string, double> doubleProperty;
    std::map<std::string, bool> boolProperty;
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
//    std::map<std::string, DoFunction> doFunction_;
    std::map<std::string, std::map<std::string, int>> settings_;

    // Objects
    std::vector<std::string> requiredObjects_;
    std::vector<Residue> residues_;
    Frame    *frame_ = nullptr;
    Frame    *cgFrame_ = nullptr;
    BondSet  *bondSet_ = nullptr;
    CGMap    *cgMap_ = nullptr;
    RDF      *rdf_ = nullptr;
    FieldMap *field_ = nullptr;
    Membrane *membrane_ = nullptr;

    // Progress updates
    const int updateFreq_[10] = {1, 2, 5, 10, 20, 50, 100, 200, 500, 1000};
    int updateLoc_ = 0;

    // Protected functions
    /** \brief Read config file and determine which functions should be performed */
    void findDoFunctions();

    /** \brief Read residues from config file - will be modified by frame later */
    void getResidues();

    /** \brief Construct objects which are required to perform requested functions */
    void setupObjects();

    /** \brief Prepare for and run the main calculation loop */
    void doMainLoop();

    /** \brief Update progress timer within the main loop */
    void updateProgress();

    /** \brief Function executed within the main loop - performs most significant work*/
    void mainLoop();

    /** \brief Perform final calculations and end program */
    void postProcess();

public:
    Common();
    ~Common();

    /** \brief Set help information */
    void setHelpStrings(const std::string &version, const std::string &header,
                        const std::string &options);

    /** \brief Verify existance of required input files */
    void collectInput(const int argc, const char *argv[],
                      const std::vector<std::string> &req_files,
                      const std::vector<std::string> &opt_files);

    /** \brief Perform all calculations */
    int run();
};

#endif //CGTOOL_COMMON_H
