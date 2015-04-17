
#include <iostream>
#include <vector>
#include <ctime>

#include <sysexits.h>

#include "frame.h"
#include "cg_map.h"
#include "itp_writer.h"
#include "parser.h"
#include "small_functions.h"
#include "file_io.h"

#ifndef NO_CMD_PARSER
#include "cmd.h"
#endif

#ifdef ELECTRIC_FIELD
#include "field_map.h"
#endif

#define PROGRESS_UPDATE_FREQ 100
#define ELECTRIC_FIELD_FREQ 100

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::clock_t;


int main(const int argc, const char *argv[]){
    clock_t start = std::clock();

    const string version_string =
            "CGTOOL v0.3.196:f2fb511fbbaa";

    const string help_header =
            "CGTOOL James Graham <J.A.Graham@soton.ac.uk> University of Southampton\n\n"
            "Performs mapping from atomistic to coarse-grained molecular dynamics\n"
            "trajectories and outputs a GROMACS ITP file containing the full mapping,\n"
            "equilibrium bond parameters and force constants.\n\n"
            "Requires GROMACS XTC and ITP files for the atomistic simulation and a\n"
            "configuration file as input.  The config file provides the mapping and\n"
            "bond parameters to be calculated as well as serveral other options.\n\n"
            "Usage:\n"
            "cgtool --dir <path to files>\n"
            "cgtool --cfg <cfg file> --xtc <xtc file> --itp <itp file>\n";
    const string help_options =
            "--cfg\tCGTOOL mapping file\tcg.cfg\t0\n"
            "--xtc\tGROMACS XTC file\tmd.xtc\t0\n"
            "--itp\tGROMACS ITP file\ttopol.top\t0\n"
            "--gro\tGROMACS GRO file\tNO DEFAULT\t0\n"
            "--dir\tDirectory containing all of the above\t./\t0\n"
            "--frm\tNumber of frames to read\t-1\t2";

    // How many threads are we using?
    int num_threads = 0;
    #pragma omp parallel reduction(+: num_threads)
    {
        num_threads = 1;
    }

    // Get input files
    split_text_output(version_string, start, num_threads);
    // If not using command line parser, replace with a simple one
    // Do this so we can compile without Boost program_options
    #ifdef NO_CMD_PARSER
    const string cfgname, xtcname, topname, groname;
    if(argc > 1 && (string(argv[1]) == "-h" || string(argv[1]) == "--help")){
        cout << help_header << endl;
        exit(EX_OK);
    }
    if(argc == 3 && string(argv[1]) == "--dir"){
        string dir = string(argv[2]);
        cfgname = dir + "/cg.cfg";
        xtcname = dir + "/md.xtc";
        topname = dir + "/topol.top";
        groname = dir + "/md.gro";
    }else if(argc == 7 && string(argv[1]) == "--cfg" && string(argv[3]) == "--xtc" && string(argv[5]) == "--itp"){
        cfgname = string(argv[2]);
        xtcname = string(argv[4]);
        topname = string(argv[6]);
        groname = "nope";
    }else{
        cout << "Wrong number of arguments provided" << endl;
        exit(EX_USAGE);
    }
    #else
    CMD cmd_parser(help_header, help_options, argc, argv);
    const string cfgname = cmd_parser.getFileArg("cfg");
    const string xtcname = cmd_parser.getFileArg("xtc");
    const string topname = cmd_parser.getFileArg("itp");
    const string groname = cmd_parser.getFileArg("gro");
    #endif

    cout << "CFG file: " << cfgname << endl;
    cout << "XTC file: " << xtcname << endl;
    cout << "ITP file: " << topname << endl;
    // GRO file is not required but helps find residues
    if(!file_exists(xtcname) || !file_exists(cfgname) || !file_exists(topname)){
        cout << "Input file does not exist" << endl;
        exit(EX_NOINPUT);
    }

    // Read number of frames from config, if not found read them all
    Parser parser(cfgname);
    vector<string> tokens;
    int num_frames_max = -1;
    if(parser.getLineFromSection("frames", tokens)) num_frames_max = stoi(tokens[0]);
    if(cmd_parser.getIntArg("frm") != 0) num_frames_max = cmd_parser.getIntArg("frm");

    int numResidues = 1;
    string resname = "";
    if(parser.getLineFromSection("residues", tokens)){
        numResidues = stoi(tokens[0]);
        resname = tokens[1];
        cout << "Mapping " << numResidues << " " << resname << " residue(s)" << endl;
    }else{
        cout << "Resname to map not found in config" << endl;
    }

    // Open files and do setup
    split_text_output("Frame setup", start, num_threads);
    Frame frame(topname, xtcname, groname, resname, numResidues);

    CGMap mapping(cfgname, resname, numResidues);
    Frame cg_frame = mapping.initFrame(frame);

    bool nomap = true;
    if(!cmd_parser.getBoolArg("nomap")){
        nomap = false;
        cg_frame.setupOutput();
    }
    BondSet bond_set(cfgname);

    #ifdef ELECTRIC_FIELD
    cout << "Doing electrostatics with " << num_threads << " thread(s)" << endl;
    FieldMap field(100, 100, 100, mapping.numBeads_);
    #endif

    // Read and process simulation frames
    split_text_output("Reading frames", start, num_threads);
    start = std::clock();
    if(num_frames_max == -1){
        cout << "Reading all frames from XTC" << endl;
    }else{
        cout << num_frames_max << " frames from XTC" << endl;
    }

    int i = 1;
    // Keep reading frames until something goes wrong (run out of frames) or hit limit
    while(frame.readNext() && (i++ < num_frames_max || num_frames_max==-1)){
        // Process each frame as we read it, frames are not retained
        #ifdef UPDATE_PROGRESS
        if(i % PROGRESS_UPDATE_FREQ == 0){
            if(num_frames_max == -1){
                cout << "Read " << i << " frames\r";
            }else{
                float tmp = i / num_frames_max;
                printf("%5.2f%%", 100.f * tmp);
                cout << " |" << string(40.f * tmp, '#') << string(40.f * (1-tmp), ' ') << "|\r";
            }
            std::flush(cout);
        }
        #endif

        if(nomap == false){
            mapping.apply(frame, cg_frame);
            cg_frame.writeToXtc();
        }

        // Calculate electric field/dipoles
        #ifdef ELECTRIC_FIELD
        if(i % ELECTRIC_FIELD_FREQ == 0) field.calculate(frame, cg_frame, mapping);
        #endif

        // Calculate bonds and store in BondStructs
        if(nomap){
            bond_set.calcBondsInternal(frame);
        }else{
            bond_set.calcBondsInternal(cg_frame);
        }
    }
    cout << endl;
    cout << "Read " << i-1 << " frames" << endl;

    // Print bitrate of XTC file input - only meaningful if we read the whole file
//    float bitrate = file_size(xtcname) * CLOCKS_PER_SEC * num_threads / ((std::clock() - start) * 1e6);
//    printf("%8.3f Mbps\n", bitrate);

    // Post processing
    split_text_output("Post processing", start, num_threads);
    if(nomap == false) cg_frame.printGRO();
    bond_set.BoltzmannInversion();

    // This bit is slow - IO limited
    #ifdef OUTPUT_CSV
    bond_set.writeCSV();
    #endif

    cout << "Printing results to ITP" << endl;
    //TODO put format choice in config file or command line option
    ITPWriter itp(resname, FileFormat::GROMACS);
    itp.printAtoms(mapping, true);
    itp.printBonds(bond_set, cmd_parser.getBoolArg("fcround"));

    // Print something so I can check results by eye
    for(int i=0; i<6 && i<bond_set.bonds_.size(); i++){
        printf("%8.4f", bond_set.bonds_[i].avg_);
    }
    if(bond_set.bonds_.size() > 6) printf("  ...");
    cout << endl;

    // Final timer
    split_text_output("Finished", start, num_threads);
    return 0;
}


