
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
#include "field_map.h"

#ifndef NO_CMD_PARSER
#include "cmd.h"
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
    clock_t very_start = std::clock();

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
            "--frames\tNumber of frames to read\t-1\t2\n"
            "--csv\tOutput bond measurements to CSV\t0\t4\n"
            "--nomap\tDon't perform cg mapping\t0\t4\n"
            "--fcround\tRound force constants\t0\t4\n"
            "--field\tCalculate electric field\t0\t4";

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
    //TODO this won't work anymore - too many things added - fix or delete
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
    if(parser.getLineFromSection("frames", tokens, 1)) num_frames_max = stoi(tokens[0]);
    if(cmd_parser.getIntArg("frames") != 0) num_frames_max = cmd_parser.getIntArg("frames");

    bool do_map = !cmd_parser.getBoolArg("nomap");
    int numResidues = 1;
    string resname = "";
    if(parser.getLineFromSection("residues", tokens, 2)){
        numResidues = stoi(tokens[0]);
        resname = tokens[1];
        if(do_map) printf("Mapping %d %s residue(s)\n", numResidues, resname.c_str());
    }else{
        cout << "Residue to map not found in config" << endl;
    }

    // Open files and do setup
    split_text_output("Frame setup", start, num_threads);
    Frame frame(topname, xtcname, groname, resname, numResidues);

    Frame cg_frame(frame);
    BondSet bond_set(cfgname);
    CGMap mapping(resname, numResidues);
    if(do_map){
        mapping.fromFile(cfgname);
        mapping.initFrame(frame, cg_frame);
        cg_frame.setupOutput();
    }

    bool do_field = cmd_parser.getBoolArg("field");
    FieldMap field(1, 1, 1, 1);
    if(do_field){
        cout << "Doing electrostatics with " << num_threads << " thread(s)" << endl;
        field.init(100, 100, 100, mapping.numBeads_);
    }

    // Read and process simulation frames
    split_text_output("Reading frames", start, num_threads);
    start = std::clock();
    if(num_frames_max == -1){
        printf("Reading all frames from XTC\n");
    }else{
        printf("Reading %6i frames from XTC\n", num_frames_max);
    }

    int i = 1;
    // Keep reading frames until something goes wrong (run out of frames) or hit limit
    while(frame.readNext() && (num_frames_max == -1 || i < num_frames_max)){
        // Process each frame as we read it, frames are not retained
        #ifdef UPDATE_PROGRESS
        if(i % PROGRESS_UPDATE_FREQ == 0){
            float time = time_since(start, num_threads);
            float fps = i / time;

            if(num_frames_max == -1){
                printf("Read %6d frames @ %d FPS\r", i, int(fps));
            }else{
                float t_remain = (num_frames_max - i) / fps;
                printf("Read %9d frames @ %d FPS %6.1fs remaining\r", i, int(fps), t_remain);
            }
            std::flush(cout);
        }
        #endif

        // Calculate bonds and store in BondStructs
        if(do_map){
            mapping.apply(frame, cg_frame);
            cg_frame.writeToXtc();
            bond_set.calcBondsInternal(cg_frame);
        }else{
            bond_set.calcBondsInternal(frame);
        }

        // Calculate electric field/dipoles
        if(do_field && i % ELECTRIC_FIELD_FREQ == 0)
                field.calculate(frame, cg_frame, mapping);

        i++;
    }

    // Print some data at the end
    cout << string(80, ' ') << "\r";
    printf("Read %9d frames", i);
    float time = time_since(start, num_threads);
    float fps = i / time;
    printf(" @ %d FPS", int(fps));
    if(num_frames_max == -1){
        // Bitrate (in MiBps) of XTC input - only meaningful if we read whole file
        float bitrate = file_size(xtcname) / (time * 1024 * 1024);
        printf("%6.1f MBps", bitrate);
    }
    printf("\n");

    // Post processing
    split_text_output("Post processing", start, num_threads);
    if(do_map){
        cg_frame.printGRO();
        bond_set.BoltzmannInversion();

        cout << "Printing results to ITP" << endl;
        //TODO put format choice in config file or command line option
        ITPWriter itp(resname, FileFormat::GROMACS);
        itp.printAtoms(mapping, true);
        itp.printBonds(bond_set, cmd_parser.getBoolArg("fcround"));
    }else{
        bond_set.calcAvgs();
    }

    // This bit is slow - IO limited
    if(cmd_parser.getBoolArg("csv")) bond_set.writeCSV();

    // Print something so I can check results by eye
    for(int i=0; i<6 && i<bond_set.bonds_.size(); i++){
        printf("%8.4f", bond_set.bonds_[i].avg_);
    }
    if(bond_set.bonds_.size() > 6) printf("  ...");
    cout << endl;

    // Final timer
    split_text_output("Finished", very_start, num_threads);
    return EX_OK;
}


