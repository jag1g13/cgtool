//
// Created by james on 14/05/15.
//

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
#include "membrane.h"

#ifndef NO_CMD_PARSER
#include "cmd.h"
#endif

#define ELECTRIC_FIELD_FREQ 100

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;


int main(const int argc, const char *argv[]){
    const double very_start = start_timer();
    double start = very_start;

    const string version_string =
            "CGTOOL v0.3.225:2aa06e3c59aa";

    const string help_header =
            "CGTOOL James Graham <J.A.Graham@soton.ac.uk> University of Southampton\n\n"
            "Usage:\n"
            "cgtool --dir <path to files>\n"
            "cgtool --cfg <cfg file> --xtc <xtc file> --itp <itp file>\n";
    const string help_options =
            "--cfg\tCGTOOL mapping file\tcg.cfg\t0\n"
            "--xtc\tGROMACS XTC file\tmd.xtc\t0\n"
            "--gro\tGROMACS GRO file\tmd.gro\t0\n"
            "--dir\tDirectory containing all of the above\t./\t0\n"
            "--frames\tNumber of frames to read\t-1\t2";

    // Allow comma separators in numbers for printf
    setlocale(LC_ALL, "");

    // ##############################################################################
    // Input Collection
    // ##############################################################################

    // Get input files
    split_text_output(version_string, start);
    CMD cmd_parser(help_header, help_options, argc, argv);
    const string cfgname = cmd_parser.getFileArg("cfg");
    const string xtcname = cmd_parser.getFileArg("xtc");
    const string groname = cmd_parser.getFileArg("gro");

    cout << "CFG file: " << cfgname << endl;
    cout << "XTC file: " << xtcname << endl;
    if(!file_exists(xtcname) || !file_exists(cfgname) || !file_exists(groname)){
        cout << "Input file does not exist" << endl;
        exit(EX_NOINPUT);
    }

    // ##############################################################################
    // System Setup
    // ##############################################################################

    // Read number of frames from config, if number not found read them all
    Parser cfg_parser(cfgname);
    vector<string> tokens;
    int num_frames_max = -1;
    if(cfg_parser.getLineFromSection("frames", tokens, 1)) num_frames_max = stoi(tokens[0]);
    if(cmd_parser.getIntArg("frames") != 0) num_frames_max = cmd_parser.getIntArg("frames");

    vector<Residue> residues;
    while(cfg_parser.getLineFromSection("residues", tokens, 1)){
        residues.emplace_back(Residue());
        Residue *res = &residues.back();
        res->resname = tokens[0];

        if(tokens.size() > 1) res->num_residues = stoi(tokens[1]);
        if(tokens.size() > 2){
            res->num_atoms = stoi(tokens[2]);
            res->calc_total();
            res->populated = true;
        }
        if(tokens.size() > 3) res->ref_atom = stoi(tokens[3]);

        printf("Mapping %'6d %5s residue(s)\n",
               residues.back().num_residues, residues.back().resname.c_str());
    }

    const int num_residues = residues.size();
    bool pop_so_far = residues[0].populated;
    residues[0].start = 0;
    for(int i=1; i<num_residues; i++){
        if(pop_so_far){
            residues[i].start = residues[i - 1].total_atoms + residues[i - 1].start;
            pop_so_far = residues[i].populated;
        }
    }

    // Open files and do setup
    split_text_output("Frame setup", start);
    string topname = "";
    Frame frame(topname, xtcname, groname, residues);
    for(Residue &res : residues) res.print();

    // Setup membrane
    Membrane mem(residues);
    mem.sortBilayer(frame);
    mem.setResolution(100);

    // Read and process simulation frames
    split_text_output("Reading frames", start);
    start = start_timer();
    if(num_frames_max == -1){
        printf("Reading all frames from XTC\n");
    }else{
        printf("Reading %'6d frames from XTC\n", num_frames_max);
    }

    // ##############################################################################
    // Main loop
    // ##############################################################################

    int i = 1;
    int progress_update_freq = 10;
    double last_update = start_timer();
    // Keep reading frames until something goes wrong (run out of frames) or hit limit
    while(frame.readNext() && (num_frames_max == -1 || i < num_frames_max)){
        // Process each frame as we read it, frames are not retained
#ifdef UPDATE_PROGRESS
        if(i % progress_update_freq == 0){
            // Set time between progress updates to nice number
            const double time_since_update = end_timer(last_update);
            if(time_since_update > 0.5f){
                if(progress_update_freq >= 10) progress_update_freq /= 10;
            }else if(time_since_update < 0.05f){
                progress_update_freq *= 10;
            }

            const double time = end_timer(start);
            const double fps = i / time;

            if(num_frames_max == -1){
                printf("Read %'9d frames @ %'d FPS\r", i, int(fps));
            }else{
                const double t_remain = (num_frames_max - i) / fps;
                printf("Read %'9d frames @ %'d FPS %6.1fs remaining\r", i, int(fps), t_remain);
            }
            std::flush(cout);

            last_update = start_timer();
        }
#endif

        // Do membrane thickness calculations - separate every 50 frames
        if(i % 50 == 0){
            mem.normalize(0);
            mem.printCSV("thickness_" + std::to_string(i));
            mem.thickness(frame, true);
        }else{
            mem.thickness(frame);
        }

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
    printf(" @ %'d FPS", int(fps));
    if(num_frames_max == -1){
        // Bitrate (in MiBps) of XTC input - only meaningful if we read whole file
        const double bitrate = file_size(xtcname) / (time * 1024 * 1024);
        printf("%6.1f MBps", bitrate);
    }
    printf("\n");

    // Post processing
    split_text_output("Post processing", start);
    // Membrane thickness simulation average
    mem.normalize(0);
    printf("Membrane thickness: %5.3f\n", mem.mean());
    mem.printCSV("thickness_running");

    // Membrane thickness in final frame
    mem.setResolution(500);
    mem.thickness(frame, true);
    mem.normalize(0);
    mem.printCSV("thickness_final");

    // Final timer
    split_text_output("Finished", very_start);
    return EX_OK;
}

