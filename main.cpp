#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <gromacs/fileio/xtcio.h>
#include <rpc/rpc.h>
#include <rpc/xdr.h>

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
using std::stoi;

int setup_frame(char*, char*, t_fileio*, int*, int*, float*, matrix, rvec**, float*);
int read_frame(char*, char*, t_fileio*, int*, int*, float*, matrix, rvec**, float*);

int main(void){
    int ok;
    char groname[18] = "test_data/npt.gro";
    char xtcname[18] = "test_data/npt.xtc";
    int natoms, step;
    float time, prec;
    rvec *x;
    matrix box;
    t_fileio* xtc;
    cout << "Starting" << endl;
    cout << "GRO file: " << groname << endl;
    //std::cin >> groname;
    cout << "XTC file: " << xtcname << endl;
    //std::cin >> xtcname;
    ok = setup_frame(groname, xtcname, xtc, &natoms, &step, &time, box, &x, &prec);
    cout << natoms << " atoms" << endl;
    cout << ok << endl;
    return !ok;
}

int setup_frame(char* groname, char* xtcname, t_fileio* xtc, int* natoms, int* step,
                float* time, matrix box, rvec** x, float* prec){
    char mode[1] = {'r'};
    char mode_write[1] = {'w'};
    char res_name[5], atom_type[3];
    string num_str, line;
    int ok = 0, num_atoms, atom_num;
    float atom_charge, atom_mass;
    ifstream gro;
    gmx_bool bOK = 0;
    xtc = open_xtc(xtcname, mode);
    ok = read_first_xtc(xtc, natoms, step, time, box, x, prec, &bOK);
    gro.open(groname);
    if(gro.is_open()){
        getline(gro, line);                     // first line of gro is a name
        //getline(gro, num_str);                  // second line is the number of atoms
        gro >> atom_num;
        cout << "Num atoms in GRO " << atom_num << endl;
        //cout << "Num from GRO " << sscanf(num_str, "%d") << endl;
        Frame frame = Frame(0, *natoms, line);
        for(int i = 0; i < *natoms; i++){       // now we can read the atoms
            Atom atom;
            gro >> res_name >> atom.atom_type >> atom.atom_num;
            //cout << (*x)[i][0] << endl;
            atom.coords = (*x)[i];
            atom.charge = 1.;
            atom.mass = 1.;
            frame.atoms.push_back(atom);
            //frame.atoms[i] = new Atom;
            //frame.atoms[i] = Atom{i, atom_type, x[i], atom_charge, atom_mass};
        }
    }
    return ok && bOK;                           // return 1 if it worked
}

//int read_frame(){
//    t_fileio *xtc_out = open_xtc(xtcfile_out, mode_write);
//    int i = 1;
//    while(ok){
//        ok_out = write_xtc(xtc_out, natoms, step, time, box, x, prec);
//        ok = read_next_xtc(xtc, natoms, &step, &time, box, *x, &prec, &bOK);
//        std::cout << "Step: " << step << std::endl;
//        i++;
//    }
//}
