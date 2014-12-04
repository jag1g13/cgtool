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

#ifndef INCLUDE_GMXFIO
#define INCLUDE_GMXFIO

#include <gromacs/fileio/xtcio.h>

#endif

#include "frame.h"
#include "cg_map.h"
#include "bondset.h"

#define DEBUG true

/*things from std that get used a lot*/
using std::ifstream;
using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::clock_t;

/*prototype functions*/
vector<float> calc_bond_lens(Frame *, vector<float>);

vector<float> calc_avg(vector<vector<float>>);

void split_text_output(const char *, const clock_t);

bool file_exists(const char *);


int main(int argc, char *argv[]){
    bool output = false;
    t_fileio *xtc, *xtc_out;
    char mode[2] = {'r', 'w'};

    /* Where does the user want us to look for input files? */
    clock_t start = std::clock();
    clock_t start_time = std::clock();
    split_text_output("Identifying files", start);
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
        strcat(xtcname, "/md.xtc");
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
    split_text_output("Frame setup", start);
    Frame frame = Frame(0, 0, "");
    xtc = open_xtc(xtcname, &mode[0]);
    if(output) xtc_out = open_xtc("out.xtc", &mode[1]);
    frame.setupFrame(groname, topname, xtc);
    Frame cg_frame = Frame(&frame);
    CGMap mapping(mapname);
    mapping.initFrame(&frame, &cg_frame);
    BondSet bond_set;
    bond_set.fromFile(bndname);

    /* Keep reading frames until something goes wrong (run out of frames) */
    split_text_output("Reading frames", start);
    start = std::clock();
    int i = 0;
    vector<vector<float>> bond_lens;
    while(frame.readNext(xtc)){
        /* Process each frame as we read it, frames are not retained */
        //cg_map(&frame, &cg_frame);
        bond_lens.push_back(bond_set.calcBondLens(&frame));
        if(i % 1000 == 0){
            cout << "Read " << i << " frames\r";
            std::flush(cout);
        }
        if(output) frame.writeToXtc(xtc_out);
        //usleep(1000);
        i++;
    }
    cout << "Read " << i << " frames" << endl;

    /* close remaining files */
    close_xtc(xtc);

    /* Post processing */
    split_text_output("Post processing", start);
    //cout << "Avg bond[0] " << calc_avg(bond_lens) << endl;
    calc_avg(bond_lens);

    /* Final timer */
    split_text_output("Finished", start);
    split_text_output("Total time", start_time);
    return 0;
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
    for(vector<vector<float>>::iterator row = bond_lens.begin(); row != bond_lens.end(); ++row){
        for(int i = 0; i < width; i++){
            sum[i] += (*row)[i];
        }
    }
    cout << "Bond lengths" << endl;
    for(int i = 0; i < width; i++){
        mean[i] = sum[i] / length;
        cout << mean[i] << endl;
    }
    return mean;
}

void split_text_output(const char *name, const clock_t start){
    clock_t now = std::clock();
    if((float) (now - start) / CLOCKS_PER_SEC > 0.1){
        cout << "--------------------" << endl;
        cout << (float) (now - start) / CLOCKS_PER_SEC << " seconds" << endl;
    }
    cout << "====================" << endl;
    cout << name << endl;
    cout << "--------------------" << endl;
}

bool file_exists(const char *name){
    struct stat buffer;
    return (stat(name, &buffer) == 0);
}

