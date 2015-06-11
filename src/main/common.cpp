//
// Created by james on 28/05/15.
//

#include "common.h"

#include <iostream>
#include <vector>

#include <sysexits.h>
#include <locale.h>

#include "small_functions.h"
#include "file_io.h"

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;

Common::Common(){
    // Allow comma separators in numbers for printf
    setlocale(LC_ALL, "");
    veryStart_ = start_timer();
    sectionStart_ = veryStart_;
}

Common::~Common(){
    if(frame_) delete frame_;
    if(cgFrame_) delete cgFrame_;
    if(bondSet_) delete bondSet_;
    if(cgMap_) delete cgMap_;
    if(rdf_) delete rdf_;
    if(field_) delete field_;
    if(membrane_) delete membrane_;
}

void Common::setHelpStrings(const std::string &version, const std::string &header,
                    const std::string &options){
    versionString_ = version;
    helpHeader_ = header;
    helpOptions_ = options;
}

void Common::collectInput(const int argc, const char *argv[],
                          const std::vector<std::string> &req_files,
                          const std::vector<std::string> &opt_files){
    split_text_output(versionString_, sectionStart_);
    CMD cmd_parser(helpHeader_, helpOptions_, argc, argv);

    // Read in files
    for(const string &f : req_files) inputFiles_[f].name = cmd_parser.getFileArg(f);
    for(const string &f : opt_files) inputFiles_[f].name = cmd_parser.getFileArg(f);

    for(auto &item : inputFiles_) item.second.exists = file_exists(item.second.name);

    for(const string &f : req_files){
        if(!inputFiles_[f].exists){
            cout << "File " << inputFiles_[f].name << " does not exist" << endl;
            exit(EX_NOINPUT);
        }
    }

    for(const auto &item : inputFiles_){
        cout << item.first << ": " << item.second.name << endl;
    }

    if(cmd_parser.getIntArg("frames") != 0) numFramesMax_ = cmd_parser.getIntArg("frames");
}

int Common::run(){
    findDoFunctions();
    getResidues();
    setupObjects();
    doMainLoop();
    postProcess();
    return EX_OK;
}

void Common::findDoFunctions(){
    Parser cfg_parser(inputFiles_["cfg"].name);

    doFunction_["map"].on =
            cfg_parser.findSection("mapping");

    doFunction_["bonds"].on =
            cfg_parser.findSection("length") || cfg_parser.findSection("angle") ||
            cfg_parser.findSection("dihedral");

    doFunction_["csv"].on =
            cfg_parser.findSection("csv");
    doFunction_["csv"].intProperty["molecules"] =
            cfg_parser.getIntKeyFromSection("csv", "molecules", 10000);

    doFunction_["rdf"].on =
            cfg_parser.findSection("rdf");
    doFunction_["rdf"].freq =
            cfg_parser.getIntKeyFromSection("rdf", "freq", 1);
    doFunction_["rdf"].doubleProperty["cutoff"] =
            cfg_parser.getDoubleKeyFromSection("rdf", "cutoff", 2.);
    doFunction_["rdf"].intProperty["resolution"] =
            cfg_parser.getIntKeyFromSection("rdf", "resolution", 100);

    doFunction_["field"].on =
            cfg_parser.findSection("field");
    doFunction_["field"].freq =
            cfg_parser.getIntKeyFromSection("field", "freq", 100);
    doFunction_["field"].intProperty["resolution"] =
            cfg_parser.getIntKeyFromSection("field", "resolution", 100);
    doFunction_["field"].intProperty["export"] =
            cfg_parser.getIntKeyFromSection("field", "export", 100);

    doFunction_["mem"].on =
            cfg_parser.findSection("membrane");
    doFunction_["mem"].freq =
            cfg_parser.getIntKeyFromSection("membrane", "calculate", 1);
    doFunction_["mem"].intProperty["export"] =
            cfg_parser.getIntKeyFromSection("membrane", "export", 100);
    doFunction_["mem"].intProperty["calculate"] =
            cfg_parser.getIntKeyFromSection("membrane", "calculate", 1);
    doFunction_["mem"].intProperty["resolution"] =
            cfg_parser.getIntKeyFromSection("membrane", "resolution", 100);
    doFunction_["mem"].intProperty["blocks"] =
            cfg_parser.getIntKeyFromSection("membrane", "blocks", 4);
    doFunction_["mem"].boolProperty["header"] =
            static_cast<bool>(cfg_parser.getIntKeyFromSection("membrane", "header", 1));
}

void Common::getResidues(){
    Parser cfg_parser(inputFiles_["cfg"].name);
    vector<string> tokens;
    if(numFramesMax_ == -1 && cfg_parser.getLineFromSection("frames", tokens, -1)){
        numFramesMax_ = stoi(tokens[0]);
    }

    while(cfg_parser.getLineFromSection("residues", tokens, 1)){
        residues_.emplace_back(Residue());
        Residue *res = &residues_.back();
        res->resname = tokens[0];

        if(tokens.size() > 1) res->num_residues = stoi(tokens[1]);
        if(tokens.size() > 2){
            res->num_atoms = stoi(tokens[2]);
            res->calc_total();
            res->populated = true;
        }
        if(tokens.size() > 3) res->ref_atom_name = tokens[3];
        if(tokens.size() > 4) throw std::runtime_error("Old input file");
    }

    const int num_residues = residues_.size();
    bool pop_so_far = residues_[0].populated;
    residues_[0].start = 0;
    for(int i=1; i<num_residues; i++){
        if(pop_so_far){
            residues_[i].start = residues_[i - 1].total_atoms + residues_[i - 1].start;
            pop_so_far = residues_[i].populated;
        }
    }
}

void Common::setupObjects(){
    // Open files and do setup
    split_text_output("Frame setup", sectionStart_);
    frame_ = new Frame(inputFiles_["xtc"].name, inputFiles_["gro"].name, residues_);
    if(inputFiles_["itp"].exists) frame_->initFromITP(inputFiles_["itp"].name);
    if(inputFiles_["fld"].exists) frame_->initFromFLD(inputFiles_["fld"].name);
    for(Residue &res : residues_) res.print();

    if(doFunction_["bonds"].on)
        bondSet_ = new BondSet(inputFiles_["cfg"].name, residues_);

    if(doFunction_["map"].on){
        cgFrame_ = new Frame(*frame_);
        cgMap_ = new CGMap(residues_);
        cgMap_->fromFile(inputFiles_["cfg"].name);
        cgMap_->initFrame(*frame_, *cgFrame_);
//        cgMap_->correctLJ();
        cgFrame_->setupOutput();
    }else{
        // If not mapping make both frames the same thing
        cgFrame_ = frame_;
    }

    if(doFunction_["rdf"].on)
        rdf_ = new RDF(residues_, doFunction_["rdf"].doubleProperty["cutoff"],
                       doFunction_["rdf"].intProperty["resolution"]);

    if(doFunction_["field"].on){
        if(doFunction_["map"].on){
            field_ = new FieldMap(doFunction_["field"].intProperty["resolution"], cgMap_->numBeads_);
        }else{
            printf("ERROR: Option 'field' requires 'mapping'\n");
            exit(EX_USAGE);
        }
    }

    if(doFunction_["mem"].on){
        membrane_ = new Membrane(residues_);
        membrane_->sortBilayer(*frame_, doFunction_["mem"].intProperty["blocks"]);
        membrane_->setResolution(doFunction_["mem"].intProperty["resolution"]);
        membrane_->header_ = doFunction_["mem"].boolProperty["header"];
    }
}

void Common::doMainLoop(){
    // Read and process simulation frames
    split_text_output("Reading frames", sectionStart_);
    sectionStart_ = start_timer();

    wholeXTCFrames_ = get_xtc_num_frames(inputFiles_["xtc"].name);
    printf("Approx %'6d frames in XTC\n", wholeXTCFrames_);

    untilEnd_ = numFramesMax_ < 0;
    if(untilEnd_){
        printf("Reading all frames from XTC\n");
    }else{
        printf("Reading %'6d frames from XTC\n", numFramesMax_);
    }

    lastUpdate_ = start_timer();

    // Process each frame as we read it, frames are not retained
    while(frame_->readNext() && (untilEnd_ || currFrame_ < numFramesMax_)){
        if(currFrame_ % updateFreq_[updateLoc_] == 0) updateProgress();
        mainLoop();
        currFrame_++;
    }

    // Print some data at the end
    cout << string(80, ' ') << "\r";
    printf("Read %'9d frames", currFrame_);
    const double time = end_timer(sectionStart_);
    const double fps = currFrame_ / time;
    printf(" @ %'d FPS", static_cast<int>(fps));
    if(numFramesMax_ == -1){
        // Bitrate (in MiBps) of XTC input - only meaningful if we read whole file
        const double bitrate = file_size(inputFiles_["xtc"].name) / (time * 1024 * 1024);
        printf("%6.1f MBps", bitrate);
    }
    printf("\n");

}

void Common::mainLoop(){
    // Calculate bonds and store in BondStructs
    if(doFunction_["map"].on){
        cgMap_->apply(*frame_, *cgFrame_);
        cgFrame_->writeToXtc();
        if(doFunction_["bonds"].on) bondSet_->calcBondsInternal(*cgFrame_);

        // Calculate electric field/dipoles
        if(doFunction_["field"].on){
            if(currFrame_ % doFunction_["field"].freq == 0){
                field_->calculate(*frame_, *cgFrame_, *cgMap_);
            }
            if(currFrame_ % doFunction_["field"].intProperty["export"] == 0){
                field_->printFieldsToFile();
            }
        }
    }else{
        if(doFunction_["bonds"].on) bondSet_->calcBondsInternal(*frame_);
    }

    if(doFunction_["rdf"].on && currFrame_ % doFunction_["rdf"].freq == 0){
        rdf_->calculateRDF(*frame_);
    }

    // Membrane thickness calculations
    if(doFunction_["mem"].on){
        if(currFrame_ % doFunction_["mem"].freq == 0){
            membrane_->thickness(*frame_);
        }
        if(doFunction_["mem"].intProperty["export"] > 0 &&
           currFrame_ % doFunction_["mem"].intProperty["export"] == 0){
            membrane_->normalize(0);
            membrane_->printCSV("thickness_" + std::to_string(currFrame_));
            membrane_->printCSVCurvature("curvature_" + std::to_string(currFrame_));
            membrane_->reset();
        }
    }
}

void Common::updateProgress(){
    // Set time between progress updates to nice number
    const double time_since_update = end_timer(lastUpdate_);
    if(time_since_update > 0.5f && updateLoc_ > 0){
        updateLoc_--;
    }else if(time_since_update < 0.01f && updateLoc_ < 10){
        updateLoc_++;
    }

    const double time = end_timer(sectionStart_);
    const double fps = currFrame_ / time;
    double t_remain = (numFramesMax_ - currFrame_) / fps;
    if(numFramesMax_ < 0) t_remain = (wholeXTCFrames_ - currFrame_) / fps;

    printf("Read %'9d frames @ %'d FPS %6.1fs remaining\r",
           currFrame_, static_cast<int>(fps), t_remain);
    std::flush(cout);

    lastUpdate_ = start_timer();
}


void Common::postProcess(){
    split_text_output("Post processing", sectionStart_);
    if(doFunction_["bonds"].on){
        if(doFunction_["map"].on){
            bondSet_->BoltzmannInversion();

            const FileFormat file_format = FileFormat::GROMACS;
            const FieldFormat field_format = FieldFormat::MARTINI;

            cout << "Printing results to ITP" << endl;
            //TODO put format choice in config file or command line option
            ITPWriter itp(residues_, file_format, field_format);
            if(cgFrame_->atomHas_.lj) itp.printAtomTypes(*cgMap_);
            itp.printAtoms(*cgMap_);
            itp.printBonds(*bondSet_);
        } else{
            bondSet_->calcAvgs();
        }

        // Write out all frame bond lengths/angles/dihedrals to file
        // This bit is slow - IO limited
        if(doFunction_["csv"].on)
            bondSet_->writeCSV(doFunction_["csv"].intProperty["molecules"]);

        // Print something so to check results by eye
        for(int j=0; j<6 && j<bondSet_->bonds_.size(); j++){
            printf("%8.4f", bondSet_->bonds_[j].avg_);
        }
        if(bondSet_->bonds_.size() > 6) printf("  ...");
        cout << endl;
    }

    if(doFunction_["map"].on) cgFrame_->printGRO();
    if(doFunction_["rdf"].on) rdf_->normalize();

    if(doFunction_["mem"].intProperty["export"] < 0){
        membrane_->normalize(0);
        membrane_->printCSV("thickness_avg");
        membrane_->printCSVCurvature("curvature_final");
    }

    // Final timer
    split_text_output("Finished", veryStart_);
}

