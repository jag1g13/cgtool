#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <stdexcept>
#include <ctime>
//#include <string>
//#include <math.h>
#include <gromacs/fileio/xtcio.h>
#include <rpc/rpc.h>
//#include <rpc/xdr.h>
#include <sys/stat.h>

#include "frame.h"

//things from std that get used a lot
using std::ifstream;
using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::clock_t;

//prototype functions
int setup_frame(char*, char*, t_fileio*, Frame*, int*, int*, float*, matrix, rvec**, float*);
int read_frame(t_fileio*, Frame*, int*, int*, float*, matrix, rvec**, float*);
vector<float> calc_bond_lens(Frame*, vector<float>);
float calc_avg(vector<float>);
void cg_map(Frame*, Frame*);
inline void split_text_output(const char *, clock_t);
inline bool file_exists(const char *);

extern int read_next_xtc(t_fileio *fio,
        int natoms, int *step, real *time,
        matrix box, rvec *x, real *prec, gmx_bool *bOK);


int main(int argc, char* argv[]){
    int ok;
    char groname[40], xtcname[40];
    int natoms, step;
    float time, prec;
    clock_t start = std::clock(), now = start;
    rvec *x;
    matrix box;
    t_fileio* xtc;
    vector<float> bond_lens;
    Frame frame = Frame(0, 0, ""), cg_frame = Frame(0, 0, "");
    char mode[1] = {'r'};
    //t_fileio *xtc_out = open_xtc("xtcout.xtc", mode_write);

    /* Where does the user want us to look for GRO and XTC files? */
    split_text_output("Identifying files", start);
    if(argc < 2) {
        cout << "Using current directory" << endl;
        strcpy(groname, "npt.gro");
        strcpy(xtcname, "md.xtc");
    }else if(argc == 2) {
        cout << "Using directory provided" << endl;
        strcpy(groname, argv[1]);
        strcat(groname, "/npt.gro");
        strcpy(xtcname, argv[1]);
        strcat(xtcname, "/md.xtc");
    }else if(argc == 3) {
        cout << "Using filenames provided" << endl;
        strcpy(groname, argv[1]);
        strcpy(xtcname, argv[2]);
    }
    if(!file_exists(groname) || !file_exists(xtcname)){
        cout << "Input file does not exist" << endl;
        throw std::runtime_error("File doesn't exist");
    }
    cout << "GRO file: " << groname << endl;
    cout << "XTC file: " << xtcname << endl;

    /* Open files and do setup */
    split_text_output("Frame setup", start);
    xtc = open_xtc(xtcname, mode);
    ok = setup_frame(groname, xtcname, xtc, &frame, &natoms, &step, &time, box, &x, &prec);

    /* Keep reading frames until something goes wrong (run out of frames) */
    split_text_output("Reading frames", start);
    start = std::clock();
    while (read_frame(xtc, &frame, &natoms, &step, &time, box, &x, &prec)) {
        /* Process each frame as we read it, frames are not retained */
        //cg_map(&frame, &cg_frame);
        bond_lens = calc_bond_lens(&frame, bond_lens);
    }
    cout << "Read " << step << " frames" << endl;

    /* close remaining files */
    close_xtc(xtc);

    /* Post processing */
    split_text_output("Post processing", start);
    cout << calc_avg(bond_lens) << endl;
    return !ok;
}

int setup_frame(char* groname, char* xtcname, t_fileio* xtc, Frame* frame,
                int* natoms, int* step, float* time,
                matrix box, rvec** x, float* prec){
    /**
    * \brief Create Frame, allocate atoms and read in data from start of XTC file
    *
    * GROMACS read_first_xtc() gets data from the XTC file about the system.
    * This function uses this data to create a Frame object to process this data
    */
    char res_name[5], line[40];
    int ok = 0, atom_num;
    //float atom_charge, atom_mass;
    ifstream gro;
    gmx_bool bOK = 0;
    Atom* atom;
    ok = read_first_xtc(xtc, natoms, step, time, box, x, prec, &bOK);
    gro.open(groname);
    if(gro.is_open()){
        gro.getline(line, 40);              // first line of gro is the run name
        cout << line << endl;
        gro >> atom_num;                    // second line is the number of atoms
        cout << "XTC num atoms:" << *natoms << endl;
        cout << "GRO num atoms:" << atom_num << endl;
        if(atom_num != *natoms){
            cout << "Number of atoms declared in XTC file "
                    "is not the same as declared in GRO file" << endl;
        }
        frame->allocate_atoms(*natoms);
        for(int i = 0; i < *natoms; i++){       // now we can read the atoms
            atom = &(frame->atoms[i]);
            gro >> res_name >> atom->atom_type >> atom->atom_num;
            memcpy(atom->coords, (*x)[i], 3*sizeof(float));
            atom->charge = 0.;
            atom->mass = 1.;
        }
    gro.close();
    //cout << "Done init" << endl;
    }else{
        cout << "GRO file cannot be opened" << endl;
    }
    return ok && bOK;                           // return 1 if it worked
}

int read_frame(t_fileio* xtc, Frame* frame,
               int* natoms, int* step, float* time,
               matrix box, rvec** x, float* prec){
    /**
    * \brief Read a frame from the XTC file into an existing Frame object
    *
    * Reads a frame into a pre-setup Frame object.
    * The same Frame object should be used for each frame to save time in allocation.
    */
    Atom* atom;
    int ok_out = 1, ok = 0, bOK = 0;
    //ok_out = write_xtc(xtc_out, *natoms, *step, *time, box, *x, *prec);
    ok = read_next_xtc(xtc, *natoms, step, time, box, *x, prec, &bOK);
    for(int i = 0; i < *natoms; i++){
        atom = &(frame->atoms[i]);                      // alias to make it tidier
        memcpy(atom->coords, (*x)[i], 3*sizeof(float)); // copy coordinates into an existing Atom
    }
    return ok_out && ok && bOK;
}

vector<float> calc_bond_lens(Frame* frame, vector<float> bond_lens){
    /**
    * \brief Calculate all relevant bond lengths in Frame
    */
    bond_lens.push_back(frame->bond_length(0, 1));
    //cout << bond_lens.back() << endl;
    return bond_lens;
}

float calc_avg(vector<float> vec){
    /**
    * \brief Calculate the average of a Vector<float>
    */
    double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
    float mean = sum / vec.size();
    return mean;
}

inline void split_text_output(const char *name, clock_t start){
    clock_t now = std::clock();

    cout << "--------------------" << endl;
    cout << (float)(now-start) / CLOCKS_PER_SEC << " seconds" << endl;
    cout << "====================" << endl;
    cout << name << endl;
    cout << "--------------------" << endl;
}

inline bool file_exists(const char *name) {
    struct stat buffer;
    return (stat (name, &buffer) == 0);
}