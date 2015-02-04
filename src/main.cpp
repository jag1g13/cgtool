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
#include "bondset.h"
#include "field_map.h"
#include "itp_writer.h"

#define UPDATE_PROGRESS true
#define PROGRESS_UPDATE_FREQ 100
#define DO_ELECTRIC_FIELD false
#define ELECTRIC_FIELD_FREQ 1000

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


int main(const int argc, const char *argv[]){
    clock_t start = std::clock();
    clock_t start_time = std::clock();

    const string version_string =
            "CGTOOL v0.1.115:fcec2412d821";

    const string help_string =
            "Requires GROMACS .xtc and .top files.\n"
            "Uses a config file to set beads and measure parameters\n\n"
            "Usage:\n"
            "cgtool\t\t\t\t; Runs using GROMACS files in the current directory\n"
            "cgtool <directory>\t\t; Runs using GROMACS files in the specified directory\n"
            "cgtool <xtc> <cfg> <top>\t; Runs using specified files - you want this one\n";

    // clang doesn't like this - it doesn't seem to handle OpenMP well?
    int num_threads = 1;
//    #pragma omp parallel
//    #pragma omp master
//    {
//        num_threads = omp_get_num_threads();
//    }
//    cout << "Running with " << num_threads << " threads" << endl;

    // get commands
//    CMD cmd_parser(help_string);
//    cmd_parser.boostParse(argc, argv);

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
    cout << "TOP file: " << topname << endl;
    cout << "CFG file: " << cfgname << endl;

    // Open files and do setup
    split_text_output("Frame setup", start, num_threads);
    Frame frame = Frame(topname, xtcname);
    CGMap mapping(cfgname);
    Frame cg_frame = mapping.initFrame(frame);
    BondSet bond_set(cfgname);
//    bond_set.fromFile(cfgname);
    FieldMap field(10, 10, 10, mapping.num_beads);

    split_text_output("Reading frames", start, num_threads);
    start = std::clock();
    vector<float> tmp;
    tmp.reserve(6);
//    vector<int> show_dipoles{0, 1, 2, 3, 4, 5};

    int i = 0;
    // Keep reading frames until something goes wrong (run out of frames)
    while(frame.readNext()){
        // Process each frame as we read it, frames are not retained
        if(i % PROGRESS_UPDATE_FREQ == 0 && UPDATE_PROGRESS){
            cout << "Read " << i << " frames\r";
            std::flush(cout);
        }
        mapping.apply(frame, cg_frame);
//        cg_frame.writeToXtc("xtcout.xtc");

        // calculate electric field/dipole
        if(i % ELECTRIC_FIELD_FREQ == 0 && DO_ELECTRIC_FIELD){
            field.setupGrid(frame);
            field.setupGridContracted(frame);
            field.calcFieldMonopolesContracted(frame);
            field.calcDipolesDirect(mapping, cg_frame, frame);
//            field.calcDipolesFit(mapping, cg_frame, frame);
            field.calcFieldDipolesContracted(cg_frame);
            field.calcTotalDipole(frame);
            field.calcSumDipole();
        }

        // calculate bonds and store in BondStructs
        bond_set.calcBondsInternal(cg_frame);

        i++;
    }
    cout << "Read " << i << " frames" << endl;

    // Post processing
    split_text_output("Post processing", start, num_threads);
    bond_set.calcAvgs();
    bond_set.writeCSV();

    ITPWriter itp("out.itp");
    itp.printAtoms(mapping);
    itp.printBonds(bond_set);

    // print something so I can tell it's working - can be removed later
    for(BondStruct bond : bond_set.bonds_){
        printf("%8.4f", bond.avg);
    }
    cout << endl;

    // Final timer
    split_text_output("Total time", start_time, num_threads);
    return 0;
}


void split_text_output(const string name, const clock_t start, const int num_threads){
    clock_t now = std::clock();
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

