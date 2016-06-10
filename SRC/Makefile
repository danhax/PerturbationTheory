
include Makefile.header

SRCS= parameters.o gaussq.o eigen.o d1mach.o utils.o getparams.o pulse.o $(ARGO) 

DEPEND= Makefile.header Makefile parameters.o

default: pert

clean:
	rm *.o *.mod *~  pert

STUFFS= pert.o $(SRCS)

pert:   $(STUFFS) 
	$(F90) $(LOADFLAGS) -o pert $(STUFFS) $(LAPACK)

parameters.o: Makefile.header Makefile parameters.F90 
	$(F90) $(FFLAGS) $(MPIFLAG) -c parameters.F90

pert.o: $(DEPEND) pert.F90 pulse.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c pert.F90 

pulse.o: $(DEPEND) pulse.F90  
	$(F90) $(FFLAGS) $(MPIFLAG) -c pulse.F90 

eigen.o: eigen.F90
	$(F90) $(FFLAGS) $(MPIFLAG) -c eigen.F90 

utils.o: $(DEPEND) utils.F90
	$(F90) $(FFLAGS) $(MPIFLAG) -c utils.F90 

getparams.o: parameters.o getparams.F90 pulse.o
	$(F90) $(FFLAGS) $(MPIFLAG) $(SAFEFLAG) $(PGFFLAG) -c getparams.F90 

####  F77 #####

d1mach.o: d1mach.f 
	$(FC) $(F77FLAGS) -c d1mach.f

gaussq.o: gaussq.f 
	$(FC) $(F77FLAGS) -c gaussq.f

arg.o: arg.c 
	$(CC) $(PGFFLAG) -c arg.c

