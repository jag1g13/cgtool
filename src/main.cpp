#include <iostream>
#include <vector>
#include <ctime>

#include <sys/stat.h>

#include "frame.h"
#include "cg_map.h"
#include "itp_writer.h"
#include "parser.h"

#ifndef NO_CMD_PARSER
#include "cmd.h"
#endif

#ifdef ELECTRIC_FIELD
#include "field_map.h"
#endif

#define PROGRESS_UPDATE_FREQ 100
#define ELECTRIC_FIELD_FREQ 100

/* things from std that get used a lot */
using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::clock_t;

/* prototype functions */
void split_text_output(const string, const clock_t, const int num_threads);
bool file_exists(const string name);


int main(const int argc, const char *argv[]){
    clock_t start = std::clock();

    const string version_string =
            "CGTOOL v0.2.167:82b4f0e76304";

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
            "--cfg\tCGTOOL mapping file\ttp.config\n"
            "--xtc\tGROMACS XTC file\tmd.xtc\n"
            "--itp\tGROMACS ITP file\ttopol.top\n"
            "--dir\tDirectory containing all of the above\t.//";

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
    string cfgname, xtcname, topname;
    if(argc > 1 && (string(argv[1]) == "-h" || string(argv[1]) == "--help")){
        cout << help_header << endl;
        exit(0);
    }
    if(argc == 3){
        if(string(argv[1]) == "--dir"){
            string dir = string(argv[2]);
            cfgname = dir + "/tp.config";
            xtcname = dir + "/md.xtc";
            topname = dir + "/topol.top";
        }
    }else if(argc == 7){
        if(string(argv[1]) == "--cfg" && string(argv[3]) == "--xtc" && string(argv[5]) == "--itp"){
            cfgname = string(argv[2]);
            xtcname = string(argv[4]);
            topname = string(argv[6]);
        }
    }else{
        cout << "Wrong number of arguments provided" << endl;
        exit(-1);
    }
    #else
    CMD cmd_parser(help_header, help_options, argc, argv);
    string cfgname = cmd_parser.getFileArg("cfg");
    string xtcname = cmd_parser.getFileArg("xtc");
    string topname = cmd_parser.getFileArg("itp");
    #endif

    cout << "Running with " << num_threads << " thread(s)" << endl;
    cout << "CFG file: " << cfgname << endl;
    cout << "XTC file: " << xtcname << endl;
    cout << "ITP file: " << topname << endl;
    if(!file_exists(xtcname) || !file_exists(cfgname) || !file_exists(topname)){
        cout << "Input file does not exist" << endl;
        exit(-1);
    }

    // Read number of frames from config, if not found read them all
    Parser parser(cfgname);
    vector<string> tokens;
    int num_frames_max = -1;
    if(parser.getLineFromSection("frames", tokens)) num_frames_max = stoi(tokens[0]);

    // Open files and do setup
    split_text_output("Frame setup", start, num_threads);
    Frame frame = Frame(topname, xtcname, cfgname);
    CGMap mapping(cfgname);
    Frame cg_frame = mapping.initFrame(frame);
    cg_frame.setupOutput("out.xtc", "out.top");
    BondSet bond_set(cfgname);

    #ifdef ELECTRIC_FIELD
    FieldMap field(10, 10, 10, mapping.numBeads_);
    #endif

    // Read and process simulation frames
    split_text_output("Reading frames", start, num_threads);
    start = std::clock();
    if(num_frames_max == -1){
        cout << "Reading all frames from XTC" << endl;
    }else{
        cout << num_frames_max << " frames from XTC" << endl;
    }

    int i = 0;
    // Keep reading frames until something goes wrong (run out of frames) or hit limit
    bool okay = true;
    while(okay && (i++ < num_frames_max || num_frames_max==-1)){
        okay = frame.readNext();
        // Process each frame as we read it, frames are not retained
        #ifdef UPDATE_PROGRESS
        if(i % PROGRESS_UPDATE_FREQ == 0){
            cout << "Read " << i << " frames\r";
            std::flush(cout);
        }
        #endif
        mapping.apply(frame, cg_frame);
        cg_frame.writeToXtc();

        // Calculate electric field/dipole
        #ifdef ELECTRIC_FIELD
        if(i % ELECTRIC_FIELD_FREQ == 0){
            field.setupGrid(frame);
            field.setupGridContracted(frame);
            field.calcFieldMonopolesContracted(frame);
            field.calcDipolesDirect(mapping, cg_frame, frame);
//            field.calcDipolesFit(mapping, cg_frame, frame);
            field.calcFieldDipolesContracted(cg_frame);
            field.calcTotalDipole(frame);
            field.calcSumDipole();
        }
        #endif

        // Calculate bonds and store in BondStructs
        bond_set.calcBondsInternal(cg_frame);
    }
    cout << "Read " << i-1 << " frames" << endl;

    // Post processing
    split_text_output("Post processing", start, num_threads);
    cg_frame.printGRO("out.gro");
    bond_set.BoltzmannInversion();

    // This bit is slow - IO limited
    #ifdef OUTPUT_CSV
    bond_set.writeCSV();
    #endif

    cout << "Printing results to ITP" << endl;
    ITPWriter itp("out.itp", frame.resname_);
    itp.printAtoms(mapping);
    itp.printBonds(bond_set);

    // Print something so I can check results by eye - can be removed later
    for(int i=0; i<6 && i<bond_set.bonds_.size(); i++){
        printf("%8.4f", bond_set.bonds_[i].avg_);
    }
    if(bond_set.bonds_.size() > 6) printf("  ...");
    cout << endl;

    // Final timer
    split_text_output("Finished", start, num_threads);
    return 0;
}


void split_text_output(const string name, const clock_t start, const int num_threads){
    clock_t now = std::clock();
    // If time has passed, how much?  Ignore small times
    if((float) (now - start) / CLOCKS_PER_SEC > 0.1){
        cout << "--------------------" << endl;
        cout << float(now - start) / (CLOCKS_PER_SEC * 1) << " seconds" << endl;
    }
    cout << "====================" << endl;
    cout << name << endl;
    cout << "--------------------" << endl;
}

bool file_exists(const string name){
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

