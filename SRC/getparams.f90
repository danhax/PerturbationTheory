
!!! MAIN SUBROUTINE FOR READING NAMELIST INPUT AND COMMAND LINE OPTIONS
  
#include "Definitions.INC"
 

subroutine getinpfile()
  use parameters
  implicit none

  integer :: nargs, getlen, i, len
#ifdef PGFFLAG
  integer :: myiargc
#endif
  character (len=SLN) :: buffer

#ifdef PGFFLAG
  nargs=myiargc()
#else
  nargs=iargc()
#endif
  do i=1,nargs
     buffer=nullbuff;     call getarg(i,buffer);     len=getlen(buffer)
     if (buffer(1:4) .eq. 'Inp=') then
        inpfile=nullbuff;        inpfile(1:len-4)=buffer(5:len)
        OFLWR "Inpfile is ", inpfile(1:len-4+1); CFL
     endif
  enddo
end subroutine getinpfile


subroutine getparams()
  use parameters
  implicit none
  integer :: nargs, getlen, i, len,  ishell, ispf,j, myiostat, iiflag,needpulse,ipulse
#ifdef PGFFLAG
  integer :: myiargc
#endif
  character (len=SLN) :: buffer
  integer :: shelltop(100)=-1            !! greater than zero: not full CI.  Number of orbitals in core shell.  Must be even.
  integer :: avectorhole(1000)=-1001
  integer :: avectorexcitefrom(1000)=-1001
  integer :: avectorexciteto(1000)=-1001
  real*8 :: tempreal,mymax

!! DUMMIES
  integer :: restrictms=0,  dfrestrictflag=0, allspinproject=1

  NAMELIST/parinp/ numstates,stateEnergies,couplingmat,expotol,par_timestep,numsteps,notiming

  OFL
  write(mpifileptr, *)
  write(mpifileptr, *) " *************************  COMMAND LINE OPTIONS  ***************************"
  write(mpifileptr, *) 
  CFL

#ifdef PGFFLAG
  nargs=myiargc()
#else
  nargs=iargc()
#endif

  open(971,file=inpfile, status="old", iostat=myiostat)
  if (myiostat/=0) then
     OFLWR "No Input.Inp found for parinp, iostat=",myiostat; CFL
  else
     OFLWR "Reading ",inpfile; CFL
     read(971,nml=parinp,iostat=myiostat)
     call checkiostat(myiostat," reading namelist PARINP")
     close(971)
  endif

  !!   NOW TAKE MCTDHF NAMELIST AND COMMAND LINE INPUT

  call openfile()

  !! NOTE THAT LBIG IS NOT A COMMAND LINE OPTION; IS IN LOOP IN GETH2OPTS.  So as of now it is set.  change later.

  do i=1,nargs
     buffer=nullbuff;     call getarg(i,buffer);     len=getlen(buffer)
     myiostat=0
     if (buffer(1:10) .eq. 'Numstates=') then
        read(buffer(11:len),*,iostat=myiostat) numstates
        write(mpifileptr,*) "Numstates set to  ", numstates, " by command line option."
     endif

     call checkiostat(myiostat," command line argument "//buffer)     
  enddo

  if (numstates.gt.MXST) then
     OFLWR "error, hardwired maximum dimension exceeded numstates MXST",numstates,MXST; CFLST
  endif

  write(mpifileptr, *) " ****************************************************************************"     
  write(mpifileptr,*);  call closefile()

  call getpulse(0)
 
end subroutine getparams


subroutine getpulse(no_error_exit_flag)   !! if flag is 0, will exit if &pulse is not read
  use parameters
  use pulse_parameters
  use pulsesubmod
  implicit none
  NAMELIST /pulse/ omega,pulsestart,pulsestrength, velflag, omega2,phaseshift,intensity,pulsetype, &
       pulsetheta,pulsephi, longstep, numpulses, chirp, ramp
  real*8 ::  time,   lastfinish, fac, pulse_end, estep
  complex*16 :: pots1(3),pots2(3),pots3(3), pots4(3), pots5(3), csumx,csumy,csumz
  integer :: i, myiostat, ipulse,no_error_exit_flag
  character (len=12) :: line
  real*8, parameter :: epsilon=1d-4
  integer, parameter :: neflux=10000
  complex*16 :: lenpot(0:neflux,3),velpot(0:neflux,3)
  real*8 :: pulseftsq(0:neflux), vpulseftsq(0:neflux)

  open(971,file=inpfile, status="old", iostat=myiostat)
  if (myiostat/=0) then
     OFLWR "No Input.Inp found, not reading pulse. iostat=",myiostat; CFL
  else
     read(971,nml=pulse,iostat=myiostat)
     if (myiostat.ne.0.and.no_error_exit_flag.eq.0) then
        OFLWR "Need &pulse namelist input!!"; CFLST
     endif
  endif
  close(971)
  OFLWR "Gauge is.... "

  select case (velflag)
  case (0)
     write(mpifileptr, *) line,"   length."
  case(1)
     write(mpifileptr, *) line,"   velocity, usual way."
  case(2)
     write(mpifileptr, *) line,"   velocity, DJH way."
  end select
  if (no_error_exit_flag.ne.0) then    !! mcscf.  just need velflag.
     return
  endif
  write(mpifileptr, *) ;  write(mpifileptr, *) "NUMBER OF PULSES:  ", numpulses;  write(mpifileptr, *) 

  lastfinish=0.d0
  do ipulse=1,numpulses

     write(mpifileptr, *) "    -----> Pulse ", ipulse," : "

     if (pulsetype(ipulse).eq.1) then
        fac=omega(ipulse)
     else
        fac=omega2(ipulse)
     endif
     if (intensity(ipulse).ne.-1.d0) then !! overrides pulsestrength
        pulsestrength(ipulse) = sqrt(intensity(ipulse)/3.51)/fac
     else
        intensity(ipulse) = (fac*pulsestrength(ipulse))**2 * 3.51  !! just output
     endif
     select case (pulsetype(ipulse))
     case (1)
        write(mpifileptr,*) "Pulse type is 1: single sine squared envelope"
        write(mpifileptr, *) "Omega, pulsestart, pulsefinish, pulsestrength:"
        write(mpifileptr, '(8F18.12)') omega(ipulse), pulsestart(ipulse), pulsestart(ipulse) + pi/omega(ipulse), pulsestrength(ipulse)
     case (2,3)
        write(mpifileptr,*) "Pulse type is 2 or 3: envelope with carrier"
        write(mpifileptr, *) "   chirp:           ", chirp(ipulse)
        write(mpifileptr, *) "   ramp:           ", ramp(ipulse), " Hartree"
        write(mpifileptr, *) "   Envelope omega:  ", omega(ipulse)
        write(mpifileptr, *) "   Pulse omega:     ", omega2(ipulse)
        write(mpifileptr, *) "   Pulsestart:      ",pulsestart(ipulse)
        write(mpifileptr, *) "   Pulsestrength:   ",pulsestrength(ipulse)
        write(mpifileptr, *) "   Intensity:       ",intensity(ipulse), " x 10^16 W cm^-2"
        write(mpifileptr, *) "   Pulsetheta:      ",pulsetheta(ipulse)
        write(mpifileptr, *) "   Pulsephi:      ",pulsephi(ipulse)
        write(mpifileptr, *) "   Pulsefinish:     ",pulsestart(ipulse) + pi/omega(ipulse)
        if (pulsetype(ipulse).eq.3) then
           write(mpifileptr,*) "---> Pulsetype 3; longstep = ", longstep(ipulse)
        else if (pulsetype(ipulse).eq.2) then
           write(mpifileptr,*) "---> Pulsetype 2."
        endif
     case (4)
        WRFL "Pulse type is 4, cw"
        write(mpifileptr, *) "   Duration omega:  ", omega(ipulse)
        write(mpifileptr, *) "   Pulse omega:     ", omega2(ipulse)
        write(mpifileptr, *) "   Pulsestrength:   ",pulsestrength(ipulse)
        write(mpifileptr, *) "   Intensity:       ",intensity(ipulse), " x 10^16 W cm^-2"
        write(mpifileptr, *) "   Pulsetheta:      ",pulsetheta(ipulse)
        write(mpifileptr, *) "   Pulsephi:      ",pulsephi(ipulse)

     end select
     if (lastfinish.lt.pulsestart(ipulse)+pi/omega(ipulse)) then
        lastfinish=pulsestart(ipulse)+pi/omega(ipulse)
     endif
  end do


end subroutine getpulse
