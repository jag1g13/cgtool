#include <iostream>
#include <vector>
#include <ctime>

#include <sys/stat.h>

#include "frame.h"
#include "cg_map.h"
#include "itp_writer.h"
#include "parser.h"

#ifdef CMD_PARSER
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
    clock_t start_time = std::clock();

    const string version_string =
            "CGTOOL v0.2.158:e23f94a3a897";

    const string help_header =
            "Requires GROMACS .xtc and .top files.\n"
            "Uses a config file to set beads and measure parameters\n\n"
            "Usage:\n";
    const string help_options =
            "--xtc\tGROMACS xtc file\tmd.xtc\n"
            "--itp\tGROMACS itp file\ttopol.top\n"
            "--cfg\tCGTOOL mapping file\tcg.cfg";
    const string help_string = help_header + help_options;

    // How many threads are we using?
    int num_threads = 0;
    #pragma omp parallel reduction(+: num_threads)
    {
        num_threads = 1;
    }

    // Get commands
    #ifdef CMD_PARSER
    CMD cmd_parser(help_options, argc, argv);
    #endif

    // Where does the user want us to look for input files?
    split_text_output(version_string, start, num_threads);
    string xtcname, topname, cfgname;
    if(argc < 2){
        cout << "Using current directory" << endl;
        xtcname = "md.xtc";
        cfgname = "tp.config";
        topname = "topol.top";
    }else if(argc == 2){
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
    }else if(argc == 4){
        cout << "Using filenames provided" << endl;
        xtcname = string(argv[1]);
        cfgname = string(argv[2]);
        topname = string(argv[3]);
    }else{
        cout << "Wrong number of arguments given" << endl;
        exit(0);
    }
    if(!file_exists(xtcname) || !file_exists(cfgname) || !file_exists(topname)){
        cout << "Input file does not exist" << endl;
        exit(-1);
    }
    cout << "Running with " << num_threads << " thread(s)" << endl;
    cout << "XTC file: " << xtcname << endl;
    cout << "CFG file: " << cfgname << endl;
    cout << "TOP file: " << topname << endl;

    // Read from config
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
    cout << num_frames_max << " frames max" << endl;
    int i = 0;
    // Keep reading frames until something goes wrong (run out of frames)
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
    }
    cout << "Read " << i-1 << " frames" << endl;

    // Post processing
    split_text_output("Post processing", start, num_threads);
    cg_frame.printGRO("out.gro");

    bond_set.calcAvgs();
    bond_set.writeCSV();
    #ifdef BOLTZMANN_INVERSION
    bond_set.boltzmannInversion();
    #endif

    cout << "Printing results to ITP" << endl;
    ITPWriter itp("out.itp");
    itp.printAtoms(mapping);
    itp.printBonds(bond_set);

    // print something so I can check results by eye - can be removed later
    for(int i=0; i<6 && i<bond_set.bonds_.size(); i++){
        printf("%8.4f", bond_set.bonds_[i].avg_);
    }
    if(bond_set.bonds_.size() > 6) printf("  ...");
    cout << endl;

    // Final timer
    split_text_output("Total time", start_time, num_threads);
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

