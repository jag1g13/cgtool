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
    // Check if frame_ and cgFrame_ point to the same thing
    if(frame_){
        if(cgFrame_ && cgFrame_ != frame_){
            delete cgFrame_;
        }
        delete frame_;
    }

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

    settings_["map"]["on"] =
            cfg_parser.findSection("mapping");

    settings_["bonds"]["on"] =
            cfg_parser.findSection("length") || cfg_parser.findSection("angle") ||
            cfg_parser.findSection("dihedral");

    settings_["csv"]["on"] =
            cfg_parser.findSection("csv");
    settings_["csv"]["molecules"] =
            cfg_parser.getIntKeyFromSection("csv", "molecules", 10000);

    settings_["rdf"]["on"] =
            cfg_parser.findSection("rdf");
    settings_["rdf"]["freq"] =
            cfg_parser.getIntKeyFromSection("rdf", "freq", 1);
    settings_["rdf"]["cutoff"] =
            cfg_parser.getIntKeyFromSection("rdf", "cutoff", 2);
    settings_["rdf"]["resolution"] =
            cfg_parser.getIntKeyFromSection("rdf", "resolution", 100);

    settings_["field"]["on"] =
            cfg_parser.findSection("field");
    settings_["field"]["freq"] =
            cfg_parser.getIntKeyFromSection("field", "freq", 100);
    settings_["field"]["resolution"] =
            cfg_parser.getIntKeyFromSection("field", "resolution", 100);
    settings_["field"]["export"] =
            cfg_parser.getIntKeyFromSection("field", "export", 100);

    settings_["mem"]["on"] =
            cfg_parser.findSection("membrane");
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
}

void Common::getResidues(){
    Parser cfg_parser(inputFiles_["cfg"].name);
    vector<string> tokens;

    while(cfg_parser.getLineFromSection("residues", tokens, 1)){
        residues_.emplace_back(Residue());
        Residue *res = &residues_.back();
        res->resname = tokens[0];

        const int size = static_cast<int>(tokens.size());

        if(size == 2) res->ref_atom_name = tokens[1];
        if(size > 2){
            printf("Reading residue from old style config file\n");
            if(size == 4) res->ref_atom_name = tokens[4];
        }
    }

    const int num_residues = static_cast<int>(residues_.size());
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

    frame_ = new Frame(inputFiles_["xtc"].name, inputFiles_["gro"].name, &residues_);
    if(inputFiles_["itp"].exists) frame_->initFromITP(inputFiles_["itp"].name);
    if(inputFiles_["fld"].exists) frame_->initFromFLD(inputFiles_["fld"].name);
    for(Residue &res : residues_) res.print();


    if(settings_["map"]["on"]){
        cgFrame_ = new Frame(*frame_, &cgResidues_);
        cgMap_ = new CGMap(&residues_, &cgResidues_);
        cgMap_->fromFile(inputFiles_["cfg"].name);
        cgMap_->initFrame(*frame_, *cgFrame_);
        cgFrame_->setupOutput();
        if(settings_["bonds"]["on"])
            bondSet_ = new BondSet(inputFiles_["cfg"].name, &cgResidues_);
    }else{
        // If not mapping make both frames point to the same thing
        cgFrame_ = frame_;
        if(settings_["bonds"]["on"])
            bondSet_ = new BondSet(inputFiles_["cfg"].name, &residues_);
    }

    if(settings_["rdf"]["on"])
        rdf_ = new RDF(&residues_, settings_["rdf"]["cutoff"]/100.,
                       settings_["rdf"]["resolution"]);

    if(settings_["field"]["on"]){
        if(settings_["map"]["on"]){
            field_ = new FieldMap(settings_["field"]["resolution"], &residues_, &cgResidues_);
        }else{
            printf("ERROR: Option 'field' requires 'mapping'\n");
            exit(EX_USAGE);
        }
    }

    if(settings_["mem"]["on"]){
        membrane_ = new Membrane(&residues_);
        membrane_->sortBilayer(*frame_, settings_["mem"]["blocks"]);
        membrane_->setResolution(settings_["mem"]["resolution"]);
        membrane_->header_ = static_cast<bool>(settings_["mem"]["header"]);
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
    if(settings_["map"]["on"]){
        cgMap_->apply(*frame_, *cgFrame_);
        cgFrame_->writeToXtc();
        if(settings_["bonds"]["on"]) bondSet_->calcBondsInternal(*cgFrame_);

        // Calculate electric field/dipoles
        if(settings_["field"]["on"]){
            if(currFrame_ % settings_["field"]["freq"] == 0){
                field_->calculate(*frame_, *cgFrame_, *cgMap_);
            }
            if(currFrame_ % settings_["field"]["export"] == 0){
                field_->printFieldsToFile();
            }
        }
    }else{
        if(settings_["bonds"]["on"]) bondSet_->calcBondsInternal(*frame_);
    }

    if(settings_["rdf"]["on"] && currFrame_ % settings_["rdf"]["freq"] == 0){
        rdf_->calculateRDF(*frame_);
    }

    // Membrane thickness calculations
    if(settings_["mem"]["on"]){
        if(currFrame_ % settings_["mem"]["freq"] == 0){
            membrane_->thickness(*frame_);
        }
        if(settings_["mem"]["export"] > 0){
            if(currFrame_ % settings_["mem"]["export"] == 0){
                membrane_->normalize(0);
                membrane_->printCSV("thickness_" + std::to_string(currFrame_));
//                membrane_->printCSVCurvature("curvature_" + std::to_string(currFrame_));
                membrane_->reset();
            }
        }
    }
}

void Common::updateProgress(){
    // Set time between progress updates to nice number
    const double time_since_update = end_timer(lastUpdate_);
    if(time_since_update > 0.5f && updateLoc_ > 0){
        updateLoc_--;
    }else if(time_since_update < 0.1f && updateLoc_ < 9){
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
    if(settings_["bonds"]["on"]){
        bondSet_->BoltzmannInversion();

        const FileFormat file_format = FileFormat::GROMACS;
        const FieldFormat field_format = FieldFormat::MARTINI;

        printf("Printing results to ITP");
        //TODO put format choice in config file or command line option
        ITPWriter itp(&residues_, file_format, field_format);
        if(cgFrame_->atomHas_.lj) itp.printAtomTypes(*cgMap_);

        if(settings_["map"]["on"]){
            itp.printAtoms(*cgMap_);
        } else{
            bondSet_->calcAvgs();
        }

        itp.printBonds(*bondSet_);

        // Write out all frame bond lengths/angles/dihedrals to file
        // This bit is slow - IO limited
        if(settings_["csv"]["on"])
            bondSet_->writeCSV(settings_["csv"]["molecules"]);

        // Print something so to check results by eye
        for(int j=0; j<6 && j<bondSet_->bonds_.size(); j++){
            printf("%8.4f", bondSet_->bonds_[j].avg_);
        }
        if(bondSet_->bonds_.size() > 6) printf("  ...");
        cout << endl;
    }

    if(settings_["map"]["on"]) cgFrame_->printGRO();
    if(settings_["rdf"]["on"]) rdf_->normalize();

    if(settings_["mem"]["export"] < 0){
        membrane_->normalize(0);
        membrane_->printCSV("thickness_avg");
        membrane_->printCSVCurvature("curvature_final");
    }

    // Final timer
    split_text_output("Finished", veryStart_);
}

