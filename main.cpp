#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
//#include <string>
//#include <math.h>
#include <gromacs/fileio/xtcio.h>
#include <rpc/rpc.h>
//#include <rpc/xdr.h>

#include "frame.h"

///* All functions return 1 if successful, 0 otherwise
// * bOK tells if a frame is not corrupted
// */
//
//t_fileio *open_xtc(const char *filename, const char *mode);
///* Open a file for xdr I/O */
//
//void close_xtc(t_fileio *fio);
///* Close the file for xdr I/O */
//
//int read_first_xtc(t_fileio *fio,
//                   int *natoms, int *step, real *time,
//                   matrix box, rvec **x, real *prec, gmx_bool *bOK);
///* Open xtc file, read xtc file first time, allocate memory for x */
//
//int read_next_xtc(t_fileio *fio,
//                  int natoms, int *step, real *time,
//                  matrix box, rvec *x, real *prec, gmx_bool *bOK);
///* Read subsequent frames */
//
//int write_xtc(t_fileio *fio,
//              int natoms, int step, real time,
//              matrix box, rvec *x, real prec);
///* Write a frame to xtc file */
//
//int xtc_check(const char *str, gmx_bool bResult, const char *file, int line);
//#define XTC_CHECK(s, b) xtc_check(s, b, __FILE__, __LINE__)
//
//void xtc_check_fat_err(const char *str, gmx_bool bResult, const char *file, int line);
//#define XTC_CHECK_FAT_ERR(s, b) xtc_check_fat_err(s, b, __FILE__, __LINE__)

//using namespace gmx;
using std::ifstream;
using std::string;
using std::cout;
using std::cin;
using std::endl;
//using std::stoi;

int setup_frame(char*, char*, t_fileio*, Frame*, int*, int*, float*, matrix, rvec**, float*);
int read_frame(t_fileio*, t_fileio*, Frame*, int*, int*, float*, matrix, rvec**, float*);
void calc_bond_lens(Frame*, std::vector<float>);
void calc_avg_bond_lens(std::vector<float>);

extern int read_next_xtc(t_fileio *fio,
        int natoms, int *step, real *time,
        matrix box, rvec *x, real *prec, gmx_bool *bOK);
/* Read subsequent frames */



int main(int argc, char* argv[]){
    int ok;
    //char groname[18] = "test_data/npt.gro";
    //char xtcname[18] = "test_data/npt.xtc";
    char groname[20], xtcname[20];
    int natoms, step;
    float time, prec;
    rvec *x;
    matrix box;
    t_fileio* xtc;
    std::vector<float> bond_lens;
    Frame frame = Frame(0, 0, "");
    char mode[1] = {'r'};
    char mode_write[1] = {'w'};
    t_fileio *xtc_out = open_xtc("xtcout.xtc", mode_write);
    if(argc < 3){
        cout << "Not enough arguments" << endl;
        return 1;
    }else {
        cout << "Starting" << endl;
        strcpy(groname, argv[1]);
        cout << "GRO file: " << groname << endl;
        strcpy(xtcname, argv[2]);
        cout << "XTC file: " << xtcname << endl;
        xtc = open_xtc(xtcname, mode);
        ok = setup_frame(groname, xtcname, xtc, &frame, &natoms, &step, &time, box, &x, &prec);
        /* Keep reading frames until something goes wrong (run out of frames) */
        while (read_frame(xtc, xtc_out, &frame, &natoms, &step, &time, box, &x, &prec)) {
            calc_bond_lens(&frame, bond_lens);
        }
        calc_avg_bond_lens(bond_lens);
        return !ok;
    }
}

int setup_frame(char* groname, char* xtcname, t_fileio* xtc, Frame* frame,
                int* natoms, int* step, float* time,
                matrix box, rvec** x, float* prec){
    char res_name[5];
    string line;
    int ok = 0, atom_num;
    //float atom_charge, atom_mass;
    ifstream gro;
    gmx_bool bOK = 0;
    Atom* atom;
    ok = read_first_xtc(xtc, natoms, step, time, box, x, prec, &bOK);
    cout << "Num atoms in XTC " << *natoms << endl;
    gro.open(groname);
    if(gro.is_open()){
        getline(gro, line);                     // first line of gro is a name
        //getline(gro, num_str);                  // second line is the number of atoms
        gro >> atom_num;
        cout << "Num atoms in GRO " << atom_num << endl;
        //cout << "Num from GRO " << sscanf(num_str, "%d") << endl;
        //Frame new_frame = Frame(0, *natoms, line);
        frame->allocate_atoms(*natoms);
        //frame = &new_frame;
        for(int i = 0; i < *natoms; i++){       // now we can read the atoms
            atom = &(frame->atoms[i]);
            gro >> res_name >> atom->atom_type >> atom->atom_num;
            memcpy(atom->coords, (*x)[i], 3*sizeof(float));
            atom->charge = 0.;
            atom->mass = 1.;
        }
    }
    cout << "Done init" << endl;
    return ok && bOK;                           // return 1 if it worked
}

int read_frame(t_fileio* xtc, t_fileio* xtc_out, Frame* frame,
               int* natoms, int* step, float* time,
               matrix box, rvec** x, float* prec){
    Atom* atom;
    int ok_out = 0, ok = 0, bOK = 0;
    //std::cout << "Step: " << *step << std::endl;
    ok_out = write_xtc(xtc_out, *natoms, *step, *time, box, *x, *prec);
    ok = read_next_xtc(xtc, *natoms, step, time, box, *x, prec, &bOK);
    //std::cout << ok_out << ok << bOK << std::endl;
    for(int i = 0; i < *natoms; i++){       // now we can read the atoms
        atom = &(frame->atoms[i]); // alias to make it tidier
        memcpy(atom->coords, (*x)[i], 3*sizeof(float));
    }
    //cout << frame->atoms[0].coords[0] << endl;
    return ok_out && ok && bOK;
}

void calc_bond_lens(Frame* frame, std::vector<float> bond_lens){
    bond_lens.push_back(frame->bond_length(0, 1));
    cout << bond_lens.back() << endl;
}

void calc_avg_bond_lens(std::vector<float> bond_lens){
    double sum = std::accumulate(bond_lens.begin(), bond_lens.end(), 0.0);
    cout << sum << endl;
    double mean = sum / bond_lens.size();
    cout << mean << endl;
}
