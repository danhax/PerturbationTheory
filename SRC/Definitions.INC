!!XXSNIPXX

#define OFL call openfile()
#define OFLWR call openfile(); write(mpifileptr,*)
#define WRFL write(mpifileptr,*)
#define CFL call closefile()
#define CFLST call closefile(); call mpistop

#define MXST 20
#define SLN 200

