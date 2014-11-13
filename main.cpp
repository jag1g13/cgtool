#include <iostream>
#include <vector>
#include <string>
#include <gromacs/fileio/xtcio.h>
#include <rpc/xdr.h>

///* All functions return 1 if successful, 0 otherwise */  
//
//extern int open_xtc(XDR *xd,char *filename,char *mode);
///* Open a file for xdr I/O */
//  
//extern void close_xtc(XDR *xd);
///* Close the file for xdr I/O */
//  
//extern int read_first_xtc(XDR *xd,char *filename,
//                          int *natoms,int *step,real *time,
//                          matrix box,rvec **x,real *prec);
///* Open xtc file, read xtc file first time, allocate memory for x */
//
//extern int read_next_xtc(XDR *xd,
//                         int *natoms,int *step,real *time,
//                         matrix box,rvec *x,real *prec);
///* Read subsequent frames */
//
//extern int write_xtc(XDR *xd,
//                     int natoms,int step,real time,
//                     matrix box,rvec *x,real prec);
///* Write a frame to xtc file */

t_fileio* setup_frame(char*, char*);

int main(void){
    t_fileio* ok;
    char grofile[20], xtcfile[20];
    std::cout << "Starting" << std::endl;
    std::cout << "GRO file: ";
    std::cin >> grofile;
    std::cout << "XTC file: ";
    std::cin >> xtcfile;
    ok = setup_frame(grofile, xtcfile);
    std::cout << ok << std::endl;
    return 0;
}

t_fileio* setup_frame(char* grofile, char* xtcfile){
    XDR xdrs;
    //xdrstdio_create(&xdrs, xtcfile, XDR_DECODE);
    t_fileio* ok = open_xtc(grofile, xtcfile);
    return ok;
}
