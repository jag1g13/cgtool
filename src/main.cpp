#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <stdexcept>
#include <ctime>
#include <map>
#include <algorithm>

#include <sys/stat.h>
//#include <omp.h>

#include "cmd.h"
#include "frame.h"
#include "cg_map.h"
#include "field_map.h"
#include "itp_writer.h"
#include "parser.h"

#define PROGRESS_UPDATE_FREQ 50
#define ELECTRIC_FIELD_FREQ 100

/* things from std that get used a lot */
using std::ifstream;
using std::ofstream;
using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::clock_t;

/* prototype functions */
vector<float> calc_avg(const vector<vector<float>> &vec);
void printToCSV(ofstream *file, const vector<float> &vec);
void split_text_output(const string, const clock_t, const int num_threads);
bool file_exists(const string name);
bool fix_PBC(const string name);


int main(const int argc, const char *argv[]){
    clock_t start = std::clock();
    clock_t start_time = std::clock();

    const string version_string =
            "CGTOOL v0.2.148:eaf90eb22cfd";

    const string help_header =
            "Requires GROMACS .xtc and .top files.\n"
            "Uses a config file to set beads and measure parameters\n\n"
            "Usage:\n";
    const string help_options =
            "--xtc\tGROMACS xtc file\tmd.xtc\n"
            "--itp\tGROMACS itp file\ttopol.top\n"
            "--cfg\tCGTOOL mapping file\tcg.cfg";
    const string help_string = help_header + help_options;

    // clang doesn't like this - it doesn't seem to do OpenMP functions?
    int num_threads = 1;
//    #pragma omp parallel
//    #pragma omp master
//    {
//        num_threads = omp_get_num_threads();
//    }
//    cout << "Running with " << num_threads << " threads" << endl;

    // get commands
//    CMD cmd_parser(help_options, argc, argv);

    // Where does the user want us to look for input files?
    split_text_output(version_string, start, num_threads);
    string xtcname, topname, cfgname;
    if(argc < 2){
        cout << "Using current directory" << endl;
        xtcname = "md.xtc";
        cfgname = "tp.config";
        topname = "topol.top";
    } else if(argc == 2){
        string arg_tmp = argv[1];
        if(arg_tmp == "-h" || arg_tmp == "--help"){
            cout << help_string << endl;
            exit(0);
        }
        cout << "Using directory provided" << endl;
        string dir = string(argv[1]);
        xtcname = dir + "/md.xtc";
        cfgname = dir + "/tp.config";
        topname = dir + "/topol.top";
    } else if(argc == 4){
        cout << "Using filenames provided" << endl;
        xtcname = string(argv[1]);
        cfgname = string(argv[2]);
        topname = string(argv[3]);
    } else{
        cout << "Wrong number of arguments given" << endl;
        throw std::runtime_error("Wrong number of arguments");
    }
    if(!file_exists(xtcname) || !file_exists(cfgname) || !file_exists(topname)){
        cout << "Input file does not exist" << endl;
        throw std::runtime_error("File doesn't exist");
    }
    cout << "XTC file: " << xtcname << endl;
    cout << "CFG file: " << cfgname << endl;
    cout << "TOP file: " << topname << endl;

    // Open files and do setup
    split_text_output("Frame setup", start, num_threads);
    Frame frame = Frame(topname, xtcname, cfgname);
    CGMap mapping(cfgname);
    Frame cg_frame = mapping.initFrame(frame);
    cg_frame.setupOutput("out.xtc", "out.top");
    BondSet bond_set(cfgname);
    FieldMap field(10, 10, 10, mapping.numBeads_);

    // Read from config
    Parser parser(cfgname);
    vector<string> tokens;
    int num_frames_max = -1;
    if(parser.getLineFromSection("frames", tokens)) num_frames_max = stoi(tokens[0]);
    cout << num_frames_max << " frames max" << endl;

    // Read and process simulation frames
    split_text_output("Reading frames", start, num_threads);
    start = std::clock();
    int i = 0;
    // Keep reading frames until something goes wrong (run out of frames)
    while(frame.readNext()){
        // Process each frame as we read it, frames are not retained
        #ifdef UPDATE_PROGRESS
        if(i % PROGRESS_UPDATE_FREQ == 0){
            cout << "Read " << i << " frames\r";
            std::flush(cout);
        }
        #endif
        mapping.apply(frame, cg_frame);
        cg_frame.writeToXtc();

        // calculate electric field/dipole
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

        // calculate bonds and store in BondStructs
        bond_set.calcBondsInternal(cg_frame);

        if(i == num_frames_max && num_frames_max != -1) break;
        i++;
    }
    cout << "Read " << i << " frames" << endl;

    // Post processing
    split_text_output("Post processing", start, num_threads);
    // output cg residue
    cg_frame.printGRO("cg.gro");

    bond_set.calcAvgs();
    bond_set.writeCSV();
    #ifdef BOLTZMANN_INVERSION
    bond_set.boltzmannInversion();
    #endif

    ITPWriter itp("out.itp");
    itp.printAtoms(mapping);
    itp.printBonds(bond_set);

    // print something so I can check results by eye - can be removed later
    for(const BondStruct &bond : bond_set.bonds_){
        printf("%8.4f", bond.avg_);
    }
    cout << endl;

    // Final timer
    split_text_output("Total time", start_time, num_threads);
    return 0;
}


void split_text_output(const string name, const clock_t start, const int num_threads){
    clock_t now = std::clock();
    // if time has passed, how much?  Ignore small times
    if((float) (now - start) / CLOCKS_PER_SEC > 0.1){
        cout << "--------------------" << endl;
        cout << float(now - start) / (CLOCKS_PER_SEC * num_threads) << " seconds" << endl;
    }
    cout << "====================" << endl;
    cout << name << endl;
    cout << "--------------------" << endl;
}

bool file_exists(const string name){
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

/** call trjconv to fix PBC problems until I work out how to do it myself */
bool fix_PBC(const string name){
    std::stringstream stream;
    stream << "trjconv -f " << name << ".xtc -o " << name << "_nojump.xtc -pbc nojump";
    return bool(system(stream.str().c_str()));
}
