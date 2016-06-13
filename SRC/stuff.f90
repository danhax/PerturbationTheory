
#include "Definitions.INC"

subroutine mpistop()
  use fileptrmod
  implicit none
  call waitawhile()
  OFLWR "STOP!";CFL
  stop
end subroutine mpistop

subroutine openfile()
end subroutine openfile

subroutine closefile()
end subroutine closefile


subroutine save1(nstates,atime,btime,tstep,vector,filename,absflag)
  use parameters
  use pulsesubmod
  implicit none
  integer,intent(in) :: nstates,atime,btime,absflag
  real*8,intent(in) :: tstep
  complex*16,intent(in) :: vector(nstates,atime:btime)
  character*(*), intent(in) :: filename
  integer :: itime

  if (absflag.eq.0) then
     open(1111, file=filename//".dat",status="unknown")
     do itime=0,numsteps
        write(1111,'(1000F18.10)') itime*par_timestep, vectdpot(itime*par_timestep,velflag), vector(:,itime)
     enddo
     close(1111)
  else
     open(1111, file=filename//"_abs.dat",status="unknown")
     do itime=0,numsteps
        write(1111,'(1000F18.10)') itime*par_timestep, vectdpot(itime*par_timestep,velflag), &
             abs(vector(:,itime)**2), SUM(abs(vector(:,itime)**2))
     enddo
     close(1111)
  endif
end subroutine save1
