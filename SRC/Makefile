
include Makefile.header

SRCS= parameters.o gaussq.o eigen.o d1mach.o utils.o getparams.o pulse.o pert.o exact.o stuff.o $(ARGO) 

DEPEND= Makefile.header Makefile parameters.o

default: pert exact

clean:
	rm *.o *.mod *~ pert exact

pert:   mainpert.o $(SRCS)
	$(F90) $(LOADFLAGS) -o pert mainpert.o $(SRCS) $(LAPACK)

exact:   mainexact.o $(SRCS)
	$(F90) $(LOADFLAGS) -o exact mainexact.o $(SRCS) $(LAPACK)

parameters.o: Makefile.header Makefile parameters.F90 
	$(F90) $(FFLAGS) $(MPIFLAG) -c parameters.F90

mainpert.o: $(DEPEND) mainpert.F90
	$(F90) $(FFLAGS) $(MPIFLAG) -c mainpert.F90 

mainexact.o: $(DEPEND) mainexact.F90
	$(F90) $(FFLAGS) $(MPIFLAG) -c mainexact.F90 

stuff.o: $(DEPEND) stuff.F90 pulse.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c stuff.F90 

pert.o: $(DEPEND) pert.F90 pulse.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c pert.F90 

exact.o: $(DEPEND) exact.F90 pulse.o
	$(F90) $(FFLAGS) $(MPIFLAG) -c exact.F90 

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

