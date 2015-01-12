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

#ifndef INCLUDE_GMXFIO
#define INCLUDE_GMXFIO

#include <gromacs/fileio/xtcio.h>

#endif

#include "frame.h"
#include "cg_map.h"
#include "bondset.h"
#include "field_map.h"

#define DEBUG true
#define UPDATE_PROGRESS true
#define PROGRESS_UPDATE_FREQ 10
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
vector<float> calc_avg(vector<vector<float>>);

void printToCSV(ofstream *file, const vector<float> *vec);

void split_text_output(const char *, const clock_t, const int num_threads);

bool file_exists(const char *);


int main(int argc, char *argv[]){
    bool output = false;
    t_fileio *xtc, *xtc_out;
    char mode[2] = {'r', 'w'};
    int num_threads = 0;
    clock_t start = std::clock();
    clock_t start_time = std::clock();

    #pragma omp parallel
    #pragma omp master
    {
        num_threads = omp_get_num_threads();
    }

    /* Where does the user want us to look for input files? */
    split_text_output("Identifying files", start, num_threads);
    cout << "Running with " << num_threads << " threads" << endl;
    char groname[40], xtcname[40], mapname[40], topname[40], bndname[40];
    if(argc < 2){
        cout << "Using current directory" << endl;
        strcpy(groname, "npt.gro");
        strcpy(xtcname, "md.xtc");
        strcpy(mapname, "sacc.map");
        strcpy(topname, "topol.top");
    } else if(argc == 2){
        cout << "Using directory provided" << endl;
        strcpy(groname, argv[1]);
        strcat(groname, "/npt.gro");
        strcpy(xtcname, argv[1]);
        strcat(xtcname, "/npt.xtc");
        strcpy(mapname, argv[1]);
        strcat(mapname, "/map.in");
        strcpy(topname, argv[1]);
        strcat(topname, "/topol.top");
        strcpy(bndname, argv[1]);
        strcat(bndname, "/bonds.in");
    } else if(argc == 5){
        cout << "Using filenames provided" << endl;
        strcpy(groname, argv[1]);
        strcpy(xtcname, argv[2]);
        strcpy(mapname, argv[3]);
        strcpy(topname, argv[4]);
    } else{
        cout << "Wrong number of arguments given" << endl;
        throw std::runtime_error("Wrong number of arguments");
    }
    if(!file_exists(groname) || !file_exists(xtcname) || !file_exists(mapname)){
        cout << "Input file does not exist" << endl;
        throw std::runtime_error("File doesn't exist");
    }
    cout << "GRO file: " << groname << endl;
    cout << "XTC file: " << xtcname << endl;
    cout << "TOP file: " << topname << endl;
    cout << "MAP file: " << mapname << endl;
    cout << "BND file: " << bndname << endl;

    /* Open files and do setup */
    split_text_output("Frame setup", start, num_threads);
    Frame frame = Frame(0, 0, "");
    xtc = open_xtc(xtcname, &mode[0]);
    if(output) xtc_out = open_xtc("out.xtc", &mode[1]);
    frame.setupFrame(groname, topname, xtc);
    Frame cg_frame = Frame(&frame);
    CGMap mapping(mapname);
    mapping.initFrame(&frame, &cg_frame);
    BondSet bond_set;
    bond_set.fromFile(bndname);
    FieldMap field(25, 25, 25, 6);
//    FieldMap field(10, 10, 10, mapping.num_beads);
//    FieldMap field(5, 5, 5, mapping.num_beads);

    /* Keep reading frames until something goes wrong (run out of frames) */
    split_text_output("Reading frames", start, num_threads);
    start = std::clock();
    int i = 0;
    vector<vector<float>> bond_lens, bond_angles, bond_dihedrals;
    vector<float> tmp;
    tmp.reserve(6);
    ofstream file_len("length.csv"), file_angle("angle.csv"), file_dih("dihedral.csv");
    while(frame.readNext(xtc)){
        /* Process each frame as we read it, frames are not retained */
        if(i % PROGRESS_UPDATE_FREQ == 0 && UPDATE_PROGRESS){
            cout << "Read " << i << " frames\r";
            std::flush(cout);
        }
        mapping.apply(&frame, &cg_frame);
        if(i % ELECTRIC_FIELD_FREQ == 0){
            field.setupGrid(&frame);
            field.setupGridContracted(&frame);
//            field.calcFieldMonopoles(&frame);
            field.calcFieldMonopolesContracted(&frame);
            field.calcDipolesDirect(&mapping, &cg_frame, &frame);
            field.calcFieldDipolesContracted(&cg_frame);
        }
//        tmp = bond_set.calcBondLens(&frame);
        tmp = bond_set.calcBondLens(&cg_frame);
        bond_lens.push_back(tmp);
        printToCSV(&file_len, &tmp);
//        tmp = bond_set.calcBondAngles(&frame);
        tmp = bond_set.calcBondAngles(&cg_frame);
        bond_angles.push_back(tmp);
        printToCSV(&file_angle, &tmp);
//        tmp = bond_set.calcBondDihedrals(&frame);
        tmp = bond_set.calcBondDihedrals(&cg_frame);
        bond_dihedrals.push_back(tmp);
        printToCSV(&file_dih, &tmp);
        if(output) frame.writeToXtc(xtc_out);
        //usleep(1000);
        i++;
    }
//    cout << "Read " << i << " frames" << endl;
//    cg_frame.printAtoms();
//    field.printDipoles();
//    field.printFields();

    /* close remaining files */
    file_len.close();
    file_angle.close();
    file_dih.close();
    close_xtc(xtc);

    /* Post processing */
//    split_text_output("Post processing", start, num_threads);
//    start = std::clock();
//    calc_avg(bond_lens);
//    calc_avg(bond_angles);
//    calc_avg(bond_dihedrals);

    /* Final timer */
//    split_text_output("Finished", start, num_threads);
    split_text_output("Total time", start_time, num_threads);
    return 0;
}

void printToCSV(ofstream *file, const vector<float> *vec){
    for(auto &item : *vec){
        *file << item << ',';
    }
    *file << endl;
}

vector<float> calc_avg(vector<vector<float>> bond_lens){
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
        cout << mean[i] << endl;
    }
    return mean;
}

void split_text_output(const char *name, const clock_t start, const int num_threads){
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

bool file_exists(const char *name){
    struct stat buffer;
    return (stat(name, &buffer) == 0);
}

