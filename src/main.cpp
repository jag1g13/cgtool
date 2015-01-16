#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <stdexcept>
#include <ctime>
#include <map>
#include <algorithm>

#include <rpc/rpc.h>
#include <sys/stat.h>
#include <omp.h>

//#ifndef INCLUDE_GMXFIO
//#define INCLUDE_GMXFIO
//#include <gromacs/fileio/xtcio.h>
//#endif

#include "frame.h"
#include "cg_map.h"
#include "bondset.h"
#include "field_map.h"

#define DEBUG true
#define UPDATE_PROGRESS true
#define PROGRESS_UPDATE_FREQ 50
#define ELECTRIC_FIELD_FREQ 50

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


int main(int argc, char *argv[]){
    clock_t start = std::clock();
    clock_t start_time = std::clock();

    // clang doesn't like this - it doesn't seem to do OpenMP
    int num_threads = 1;
//    #pragma omp parallel
//    #pragma omp master
//    {
//        num_threads = omp_get_num_threads();
//    }
//    cout << "Running with " << num_threads << " threads" << endl;

    // Where does the user want us to look for input files?
    split_text_output("Identifying files", start, num_threads);
    string groname, xtcname, topname, cfgname;
    if(argc < 2){
        cout << "Using current directory" << endl;
        groname = "npt.gro";
        xtcname = "md.xtc";
        cfgname = "tp.config";
        topname = "topol.top";
    } else if(argc == 2){
        cout << "Using directory provided" << endl;
        string dir = string(argv[1]);
        groname = dir + "/npt.gro";
        xtcname = dir + "/md.xtc";
        cfgname = dir + "/tp.config";
        topname = dir + "/topol.top";
    } else if(argc == 5){
        cout << "Using filenames provided" << endl;
        groname = string(argv[1]);
        xtcname = string(argv[2]);
        cfgname = string(argv[3]);
        topname = string(argv[4]);
    } else{
        cout << "Wrong number of arguments given" << endl;
        throw std::runtime_error("Wrong number of arguments");
    }
    if(!file_exists(groname) || !file_exists(xtcname) || !file_exists(cfgname)){
        cout << "Input file does not exist" << endl;
        throw std::runtime_error("File doesn't exist");
    }
    cout << "GRO file: " << groname << endl;
    cout << "XTC file: " << xtcname << endl;
    cout << "TOP file: " << topname << endl;
    cout << "CFG file: " << cfgname << endl;

    // Open files and do setup
    split_text_output("Frame setup", start, num_threads);
    Frame frame = Frame(groname, topname, cfgname, xtcname);
    CGMap mapping(cfgname);
    Frame cg_frame = mapping.initFrame(frame);
    BondSet bond_set;
    bond_set.fromFile(cfgname);
    FieldMap field(10, 10, 10, mapping.num_beads);

    split_text_output("Reading frames", start, num_threads);
    start = std::clock();
    vector<vector<float>> bond_lens, bond_angles, bond_dihedrals;
    vector<float> tmp;
    tmp.reserve(6);
//    vector<int> show_dipoles{0, 1, 2, 3, 4, 5};
    ofstream file_len("length.csv"), file_angle("angle.csv"), file_dih("dihedral.csv");
    int i = 0;
    // Keep reading frames until something goes wrong (run out of frames)
    while(frame.readNext()){
        // Process each frame as we read it, frames are not retained
        if(i % PROGRESS_UPDATE_FREQ == 0 && UPDATE_PROGRESS){
            cout << "Read " << i << " frames\r";
            std::flush(cout);
        }
        mapping.apply(frame, cg_frame);
        // calculate electric field/dipole
        if(i % ELECTRIC_FIELD_FREQ == 0){
            field.setupGrid(&frame);
            field.setupGridContracted(&frame);
//            field.calcFieldMonopoles(&frame);
            field.calcFieldMonopolesContracted(&frame);
            field.calcDipolesDirect(&mapping, &cg_frame, &frame);
//            field.calcDipolesFit(&mapping, &cg_frame, &frame);
            field.calcFieldDipolesContracted(&cg_frame);
//            field.calcTotalDipole(&frame);
//            field.calcSumDipole(show_dipoles);
        }

        // calculate bonds
        tmp = bond_set.calcBondLens(&cg_frame);
        bond_lens.push_back(tmp);
        printToCSV(&file_len, tmp);
        tmp = bond_set.calcBondAngles(&cg_frame);
        bond_angles.push_back(tmp);
        printToCSV(&file_angle, tmp);
        tmp = bond_set.calcBondDihedrals(&cg_frame);
        bond_dihedrals.push_back(tmp);

        printToCSV(&file_dih, tmp);
        i++;
    }
    cout << "Read " << i << " frames" << endl;
//    cg_frame.printAtoms();
//    field.printDipoles();
//    field.printFields();

    // close remaining files
    file_len.close();
    file_angle.close();
    file_dih.close();

    // Post processing
    split_text_output("Post processing", start, num_threads);
    start = std::clock();
    calc_avg(bond_lens);
    calc_avg(bond_angles);
    calc_avg(bond_dihedrals);

    // Final timer
    split_text_output("Total time", start_time, num_threads);
    return 0;
}

void printToCSV(ofstream *file, const vector<float> &vec){
    for(const auto &item : vec){
        *file << item << ',';
    }
    *file << endl;
}

vector<float> calc_avg(const vector<vector<float>> &bond_lens){
    /**
    * \brief Calculate the average of bond lengths in vector<vector<float>
    *
    * This is not cache friendly
    */
    int length = bond_lens.size();
    int width = bond_lens[0].size();
    vector<float> sum(width);
    vector<float> mean(width);
    for(auto &row : bond_lens){
        for(int i = 0; i < width; i++){
            sum[i] += row[i];
        }
    }
    cout << "Bonds" << endl;
    for(int i = 0; i < width; i++){
        mean[i] = sum[i] / length;
//        cout << mean[i] << "\t";
        printf("%8.3f", mean[i]);
    }
    cout << endl;
    return mean;
}

void split_text_output(const string name, const clock_t start, const int num_threads){
    clock_t now = std::clock();
    if((float) (now - start) / CLOCKS_PER_SEC > 0.1){
        cout << "--------------------" << endl;
//        cout << (float) (now - start) / (CLOCKS_PER_SEC) << " seconds" << endl;
        cout << (float) (now - start) / (CLOCKS_PER_SEC * num_threads) << " seconds" << endl;
    }
    cout << "====================" << endl;
    cout << name << endl;
    cout << "--------------------" << endl;
}

bool file_exists(const string name){
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

