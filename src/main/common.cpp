//
// Created by james on 28/05/15.
//

#include "common.h"

#include <iostream>
#include <vector>

#include <sysexits.h>
#include <locale.h>

#include "small_functions.h"

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

    if(cgMap_) delete cgMap_;
}

void Common::setHelpStrings(const std::string &version, const std::string &header,
                            const std::string &options, const string &compile){
    versionString_ = version;
    helpHeader_ = header;
    helpOptions_ = options;
    compileInfo_ = compile;
}

void Common::collectInput(const int argc, const char *argv[],
                          const std::vector<std::string> &req_files,
                          const std::vector<std::string> &opt_files){
    split_text_output(versionString_, sectionStart_);
    CMD cmd_parser(helpHeader_, helpOptions_, compileInfo_, argc, argv);

    // Read in files
    for(const string &f : req_files) inputFiles_[f].name = cmd_parser.getStringArg(f);
    for(const string &f : opt_files) inputFiles_[f].name = cmd_parser.getStringArg(f);

    for(auto &item : inputFiles_) item.second.exists = file_exists(item.second.name);

    for(const string &f : req_files){
        if(!inputFiles_[f].exists){
            printf("ERROR: File %s does not exist\n", inputFiles_[f].name.c_str());
            exit(EX_NOINPUT);
        }
    }

    if(cmd_parser.getIntArg("frames") != 0) numFramesMax_ = cmd_parser.getIntArg("frames");
}

int Common::run(){
    readConfig();
    getResidues();

    // Open files and do setup
    split_text_output("Frame setup", sectionStart_);
    setupObjects();

    doMainLoop();

    split_text_output("Post processing", sectionStart_);
    postProcess();

    // Final timer
    split_text_output("Finished", veryStart_);
    return EX_OK;
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
            printf("NOTE: Reading residue from old style config file\n");
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


void Common::doMainLoop(){
    // Read and process simulation frames
    split_text_output("Reading frames", sectionStart_);
    sectionStart_ = start_timer();

    wholeXTCFrames_ = get_xtc_num_frames(inputFiles_["xtc"].name);
    printf("Approx %'8d frames in XTC\n", wholeXTCFrames_);

    untilEnd_ = numFramesMax_ < 0;
    if(untilEnd_){
        printf("Reading all frames from XTC\n");
    }else{
        printf("Reading %'7d frames from XTC\n", numFramesMax_);
    }

    lastUpdate_ = start_timer();

    // Process each frame as we read it, frames are not retained
    bool end = false;
    while(!end){
        end = !(frame_->readNext() && (untilEnd_ || currFrame_ < numFramesMax_));
        if(currFrame_ % updateFreq_[updateLoc_] == 0) updateProgress();
        currFrame_++;
        mainLoop();
    }

    // Print some data at the end
    cout << string(80, ' ') << "\r";
    printf("Read %'10d frames", currFrame_);
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

void Common::updateProgress(){
    // Set time between progress updates to nice number
    const double time_since_update = end_timer(lastUpdate_);
    if(time_since_update > 0.5f && updateLoc_ > 0){
        updateLoc_--;
    }else if(time_since_update < 0.1f && updateLoc_ < 10){
        updateLoc_++;
    }

    const double time = end_timer(sectionStart_);
    const double fps = currFrame_ / time;
    double t_remain = (numFramesMax_ - currFrame_) / fps;
    if(numFramesMax_ < 0) t_remain = (wholeXTCFrames_ - currFrame_) / fps;

    printf("Read %'10d frames @ %'d FPS %6.1fs remaining\r",
           currFrame_, static_cast<int>(fps), t_remain);
    std::flush(cout);

    lastUpdate_ = start_timer();
}
