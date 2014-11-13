
/* All functions return 1 if successful, 0 otherwise */  

extern int open_xtc(XDR *xd,char *filename,char *mode);
/* Open a file for xdr I/O */
  
extern void close_xtc(XDR *xd);
/* Close the file for xdr I/O */
  
extern int read_first_xtc(XDR *xd,char *filename,
                          int *natoms,int *step,real *time,
                          matrix box,rvec **x,real *prec);
/* Open xtc file, read xtc file first time, allocate memory for x */

extern int read_next_xtc(XDR *xd,
                         int *natoms,int *step,real *time,
                         matrix box,rvec *x,real *prec);
/* Read subsequent frames */

extern int write_xtc(XDR *xd,
                     int natoms,int step,real time,
                     matrix box,rvec *x,real prec);
/* Write a frame to xtc file */