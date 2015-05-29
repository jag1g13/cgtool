//
// Created by james on 28/05/15.
//

#include "common.h"

#include <iostream>
#include <vector>

#include <sysexits.h>
#include <locale.h>

#include "frame.h"
#include "cg_map.h"
#include "itp_writer.h"
#include "parser.h"
#include "small_functions.h"
#include "file_io.h"
#include "field_map.h"
#include "cmd.h"
#include "rdf.h"

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
        cout << std::toupper(item.first) << ": " << item.second.name << endl;
    }

    if(cmd_parser.getIntArg("frames") != 0) numFramesMax_ = cmd_parser.getIntArg("frames");
}

void Common::findDoFunctions(){
    Parser cfg_parser(inputFiles_["cfg"].name);

    doFunction_["map"].on = cfg_parser.findSection("mapping");

    doFunction_["rdf"].on = cfg_parser.findSection("rdf");
    doFunction_["rdf"].freq = cfg_parser.getIntKeyFromSection("rdf", "freq", 1);

    doFunction_["field"].on = cfg_parser.findSection("field");
    doFunction_["field"].freq = cfg_parser.getIntKeyFromSection("field", "freq", 100);
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

int Common::do_stuff(const int argc, const char *argv[]){
    double start = veryStart_;

    // ##############################################################################
    // System Setup
    // ##############################################################################

    // Read number of frames from config, if number not found read them all
    int num_frames_max = -1;


    // Open files and do setup
    split_text_output("Frame setup", start);
    Frame frame(inputFiles_["xtc"].name, inputFiles_["gro"].name, residues_);
    if(inputFiles_["itp"].exists) frame.initFromITP(inputFiles_["itp"].name);
    if(inputFiles_["fld"].exists) frame.initFromFLD(inputFiles_["fld"].name);
    for(Residue &res : residues_) res.print();

    Frame cg_frame(frame);
    BondSet bond_set(inputFiles_["cfg"].name, residues_);
    CGMap mapping(residues_);
    if(doFunction_["map"].on){
        mapping.fromFile(inputFiles_["cfg"].name);
        mapping.initFrame(frame, cg_frame);
        mapping.correctLJ();
        cg_frame.setupOutput();
    }

    RDF rdf;
    if(doFunction_["rdf"].on){
        const double cutoff = cfg_parser.getDoubleKeyFromSection("rdf", "cutoff", 2.);
        const int resolution = cfg_parser.getIntKeyFromSection("rdf", "resolution", 100.);
        rdf.init(residues_, cutoff, resolution);
    }

    FieldMap field;
    if(doFunction_["field"].on) field.init(100, 100, 100, mapping.numBeads_);

    // Read and process simulation frames
    split_text_output("Reading frames", start);
    start = start_timer();
    const int full_xtc_frames = get_xtc_num_frames(inputFiles_["xtc"].name);
    printf("Total of %'6d frames in XTC\n", full_xtc_frames);
    if(num_frames_max < 0){
        printf("Reading all frames from XTC\n");
    }else{
        printf("Reading %'6d frames from XTC\n", num_frames_max);
    }

    // ##############################################################################
    // Main loop
    // ##############################################################################

    int i = 1;
    int progress_update_loc = 0;
    const int progress_update_freqs[6] = {1, 2, 5, 10, 100, 1000};
    double last_update = start_timer();
    // Keep reading frames until something goes wrong (run out of frames) or hit limit
    while(frame.readNext() && (num_frames_max < 0 || i < num_frames_max)){
        // Process each frame as we read it, frames are not retained
#ifdef UPDATE_PROGRESS
        if(i % progress_update_freqs[progress_update_loc] == 0){
            // Set time between progress updates to nice number
            const double time_since_update = end_timer(last_update);
            if(time_since_update > 0.5f && progress_update_loc > 0){
                progress_update_loc--;
            }else if(time_since_update < 0.01f && progress_update_loc < 6){
                progress_update_loc++;
            }

            const double time = end_timer(start);
            const double fps = i / time;

            double t_remain = (num_frames_max - i) / fps;
            if(num_frames_max < 0) t_remain = (full_xtc_frames - i) / fps;
            printf("Read %'9d frames @ %'d FPS %6.1fs remaining\r", i, static_cast<int>(fps), t_remain);
            std::flush(cout);

            last_update = start_timer();
        }
#endif

        // Calculate bonds and store in BondStructs
        if(doFunction_["map"].on){
            mapping.apply(frame, cg_frame);
            cg_frame.writeToXtc();
            bond_set.calcBondsInternal(cg_frame);
        }else{
            bond_set.calcBondsInternal(frame);
        }

        // Calculate electric field/dipoles
        if(doFunction_["field"].on && i % doFunction_["field"].freq == 0){
            field.calculate(frame, cg_frame, mapping);
        }

        if(doFunction_["rdf"].on && i % doFunction_["rdf"].freq == 0) rdf.calculateRDF(frame);

        i++;
    }

    // ##############################################################################
    // Post processing / Averaging
    // ##############################################################################

    // Print some data at the end
    cout << string(80, ' ') << "\r";
    printf("Read %'9d frames", i);
    const double time = end_timer(start);
    const double fps = i / time;
    printf(" @ %'d FPS", static_cast<int>(fps));
    if(num_frames_max == -1){
        // Bitrate (in MiBps) of XTC input - only meaningful if we read whole file
        const double bitrate = file_size(inputFiles_["xtc"].name) / (time * 1024 * 1024);
        printf("%6.1f MBps", bitrate);
    }
    printf("\n");

    // Post processing
    split_text_output("Post processing", start);
    if(doFunction_["map"].on){
        cg_frame.printGRO();
        bond_set.BoltzmannInversion();

        const FileFormat file_format = FileFormat::GROMACS;
        const FieldFormat field_format = FieldFormat::MARTINI;

        cout << "Printing results to ITP" << endl;
        //TODO put format choice in config file or command line option
        ITPWriter itp(residues_, file_format, field_format);
        if(cg_frame.atomHas_.lj) itp.printAtomTypes(mapping);
        itp.printAtoms(mapping);
        itp.printBonds(bond_set, cmd_parser.getBoolArg("fcround"));
    }else{
        bond_set.calcAvgs();
    }

    // Write out all frame bond lengths/angles/dihedrals to file
    // This bit is slow - IO limited
    if(cmd_parser.getBoolArg("csv")) bond_set.writeCSV();

    // Calculate RDF
    if(doFunction_["rdf"].on) rdf.normalize();

    // Print something so I can check results by eye
    for(int j=0; j<6 && j<bond_set.bonds_.size(); j++){
        printf("%8.4f", bond_set.bonds_[j].avg_);
    }
    if(bond_set.bonds_.size() > 6) printf("  ...");
    cout << endl;

    // Final timer
    split_text_output("Finished", veryStart_);
    return EX_OK;
}

