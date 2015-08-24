#include <string>
#include <vector>

#include "ramsi.h"

using std::string;
using std::vector;

int main(const int argc, const char *argv[]){
    const string version_string =
            "RAMSi v0.3"
            #include "revision_number.inc"
    ;

    const string help_header =
            "James Graham <J.A.Graham@soton.ac.uk> University of Southampton\n\n"
            "Performs several analysis functions for biomembrane simulations in GROMACS.\n\n"
            "Requires GROMACS XTC and GRO files from the simulation and a configuration\n"
            "file as input.  The config file contains the output settings and a\n"
            "coarse-grain mapping if a AA->CG transformation is required.\n\n"
            "Usage:\n"
            "ramsi -c <CFG file> -x <XTC file> -g <GRO file>\n";
    // Option syntax is <long flag> \t <comment> \t <flag type> [ \t <default value>]
    // Flag types are 0 - path, 1 - string, 2 - int, 3 - float, 4 - bool
    const string help_options =
            "--cfg\tRAMSi config file\t0\n"
            "--xtc\tGROMACS XTC file\t0\n"
            "--gro\tGROMACS GRO file\t0";

    const string compile_info =
            #include "compile_info.inc"
    ;

    Ramsi ramsi;
    ramsi.setHelpStrings(version_string, help_header, help_options, compile_info);

    vector<string> req_files = {"cfg", "xtc", "gro"};
    vector<string> opt_files = {};

    ramsi.collectInput(argc, argv, req_files, opt_files);
    return ramsi.run();
}


void Ramsi::readConfig(){
    Parser cfg_parser(inputFiles_["cfg"].name);

    settings_["map"]["on"] =
            cfg_parser.findSection("mapping");

    settings_["mem"]["freq"] =
            cfg_parser.getIntKeyFromSection("membrane", "calculate", 1);
    settings_["mem"]["export"] =
            cfg_parser.getIntKeyFromSection("membrane", "export", 100);
    settings_["mem"]["calculate"] =
            cfg_parser.getIntKeyFromSection("membrane", "calculate", 1);
    settings_["mem"]["resolution"] =
            cfg_parser.getIntKeyFromSection("membrane", "resolution", 100);
    settings_["mem"]["blocks"] =
            cfg_parser.getIntKeyFromSection("membrane", "blocks", 4);
    settings_["mem"]["header"] =
            cfg_parser.getIntKeyFromSection("membrane", "header", 1);

    numFramesMax_ = cfg_parser.getIntKeyFromSection("general", "frames", -1);
}

void Ramsi::setupObjects(){
    frame_ = new Frame(inputFiles_["xtc"].name, inputFiles_["gro"].name, &residues_);
    for(Residue &res : residues_) res.print();

    if(settings_["map"]["on"]){
        cgFrame_ = new Frame(*frame_, &cgResidues_);
        cgMap_ = new CGMap(&residues_, &cgResidues_);
        cgMap_->fromFile(inputFiles_["cfg"].name);
        cgMap_->initFrame(*frame_, *cgFrame_);
//        cgFrame_->setupOutput();
//        cgFrame_->outputSingleFrame();
    }else{
        // If not mapping make both frames point to the same thing
        cgFrame_ = frame_;
    }

    membrane_ = new Membrane(&residues_);
    membrane_->sortBilayer(*frame_, settings_["mem"]["blocks"]);
    membrane_->setResolution(settings_["mem"]["resolution"]);
    membrane_->header_ = static_cast<bool>(settings_["mem"]["header"]);
}

void Ramsi::mainLoop(){
    if(settings_["map"]["on"]){
        cgMap_->apply(*frame_, *cgFrame_);
//        cgFrame_->outputTrajectoryFrame();
    }

    // Membrane calculations
    if(currFrame_ % settings_["mem"]["freq"] == 0){
        membrane_->thickness(*cgFrame_);
        membrane_->curvature(*cgFrame_);
    }

    if(settings_["mem"]["export"] > 0){
        if(currFrame_ % settings_["mem"]["export"] == 0){
            membrane_->normalize(0);
            membrane_->printCSV("thickness_" + std::to_string(currFrame_));
            membrane_->printCSVCurvature("curvature_" + std::to_string(currFrame_));
            membrane_->reset();
        }
    }
}

void Ramsi::postProcess(){
    if(settings_["mem"]["export"] < 0){
        membrane_->normalize(0);
        membrane_->printCSV("thickness_avg");
        membrane_->printCSVCurvature("curvature_final");
    }
}

Ramsi::~Ramsi(){
    if(membrane_) delete membrane_;
}
