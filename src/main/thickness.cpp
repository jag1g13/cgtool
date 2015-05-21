//
// Created by james on 14/05/15.
//

#include <iostream>
#include <vector>

#include <sysexits.h>
#include <locale.h>

#include "frame.h"
#include "parser.h"
#include "small_functions.h"
#include "membrane.h"
#include "cmd.h"

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::stoi;


int main(const int argc, const char *argv[]){
    const double very_start = start_timer();
    double start = very_start;

    const string version_string =
            "CGTOOL THICKNESS v0.3.232:63326f510f40";

    const string help_header =
            "CGTOOL THICKNESS James Graham <J.A.Graham@soton.ac.uk> University of Southampton\n\n"
            "Usage:\n"
            "cgtool_thickness -c <cfg file> -x <xtc file> -g <gro file> [-f <num frames>]\n";
    const string help_options =
            "--cfg\tCGTOOL config file\t0\n"
            "--xtc\tGROMACS XTC file\t0\n"
            "--gro\tGROMACS GRO file\t0\n"
            "--dir\tDirectory containing all of the above\t0\n"
            "--frames\tNumber of frames to read\t2\t-1\n"
            "--header\tPrint file header in membrane export\t4\t0";

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
    cout << "GRO file: " << groname << endl;
    if(!file_exists(cfgname)){
        cout << "CFG file does not exist" << endl;
        exit(EX_NOINPUT);
    }
    if(!file_exists(xtcname)){
        cout << "XTC file does not exist" << endl;
        exit(EX_NOINPUT);
    }
    if(!file_exists(groname)){
        cout << "GRO file does not exist" << endl;
        exit(EX_NOINPUT);
    }


    // ##############################################################################
    // System Setup
    // ##############################################################################

    // Read number of frames from config or command line, if number not found read them all
    Parser cfg_parser(cfgname);
    vector<string> tokens;
    int num_frames_max = -1;
    if(cfg_parser.getLineFromSection("frames", tokens, 1)) num_frames_max = stoi(tokens[0]);
    if(cmd_parser.getIntArg("frames") != 0) num_frames_max = cmd_parser.getIntArg("frames");

    // Get residues from config
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
    Frame frame(xtcname, groname, residues);

    // Print residue info
    printf("\nResidues\n----------\n");
    for(Residue &res : residues) res.print();

    // Read settings from config or default value
    const int exp_every_N  = cfg_parser.getIntKeyFromSection("membrane", "export",     100);
    const int calc_every_N = cfg_parser.getIntKeyFromSection("membrane", "calculate",    1);
    const int resolution   = cfg_parser.getIntKeyFromSection("membrane", "resolution", 100);
    const int blocks       = cfg_parser.getIntKeyFromSection("membrane", "blocks",       4);

    // Setup membrane
    Membrane mem(residues);
    printf("\nBilayer\n----------\n");
    mem.sortBilayer(frame, blocks);
    mem.setResolution(resolution);
    const bool mem_header = cmd_parser.getBoolArg("header");
    if(mem_header) printf("Exporting membrane thickness with header\n");

    // Calculate from GRO file
    mem.thickness(frame);
    mem.normalize(0);
    mem.printCSV("thickness_gro", mem_header);
    mem.reset();

    // Read and process simulation frames
    split_text_output("Reading frames", start);
    start = start_timer();
    const int full_xtc_frames = get_xtc_num_frames(xtcname);
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
                if(progress_update_freq >= 2) progress_update_freq /= 2;
            }else if(time_since_update < 0.2f){
                progress_update_freq *= 2;
            }

            const double time = end_timer(start);
            const double fps = i / time;

            double t_remain = (num_frames_max - i) / fps;
            if(num_frames_max < 0) t_remain = (full_xtc_frames - i) / fps;
            printf("Read %'9d frames @ %'d FPS %6.1fs remaining\r", i, int(fps), t_remain);
            std::flush(cout);

            last_update = start_timer();
        }
#endif

        // Do membrane thickness calculations - separate every N frames
        if(i % calc_every_N == 0) mem.thickness(frame);
        if(exp_every_N > 0 && i % exp_every_N == 0){
            mem.normalize(0);
            mem.printCSV("thickness_" + std::to_string(i), mem_header);
            mem.reset();
        }

        i++;
    }

    // ##############################################################################
    // Post processing / Averaging
    // ##############################################################################

    // Print timing data at the end
    cout << string(80, ' ') << "\r";
    const double time = end_timer(start);
    const double MBrate = file_size(xtcname) / (time * 1024 * 1024);
    printf("Read %'9d frames @ %'4d FPS %'6.1f MBps\n", i, int(i/time), MBrate);

    // Calculate thickness average if requested
    if(exp_every_N < 0){
        mem.normalize(0);
        mem.printCSV("thickness_average", mem_header);
    }

    // Final timer
    split_text_output("Finished", very_start);
    return EX_OK;
}

