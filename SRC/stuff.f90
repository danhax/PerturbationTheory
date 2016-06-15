
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


subroutine save1(atime,btime,tstep,vector,filename,absflag)
  use parameters
  use pulsesubmod
  implicit none
  integer,intent(in) :: atime,btime,absflag
  real*8,intent(in) :: tstep
  complex*16,intent(in) :: vector(numstates,atime:btime)
  complex*16 :: facs(numstates)
  character*(*), intent(in) :: filename
  integer :: itime

  if (absflag.eq.0) then
     open(1111, file=filename//".dat",status="unknown")
     do itime=atime,btime
        facs(:)=exp((0d0,+1d0)*itime*par_timestep*myStateEnergies(1:numstates))
        write(1111,'(1000F18.10)') itime*par_timestep, vectdpot(itime*par_timestep,velflag), &
             vector(1:min(numstates,maxoutstates),itime)*facs(1:min(numstates,maxoutstates))
     enddo
     close(1111)
  else
     open(1111, file=filename//"_abs.dat",status="unknown")
     do itime=atime,btime
        write(1111,'(1000F18.10)') itime*par_timestep, vectdpot(itime*par_timestep,velflag), &
             abs(vector(1:min(numstates,maxoutstates),itime)**2), SUM(abs(vector(:,itime)**2))
     enddo
     close(1111)
  endif
end subroutine save1


subroutine get_velocityop(avector,velocityop)
  use parameters
  use expsubmod
  implicit none
  real*8,intent(in) :: avector !! vector potential
  complex*16,intent(out) :: velocityop(numstates,numstates)
  complex*16 :: gaugetransform(numstates,numstates),&
       h0mat(numstates,numstates)   !! AUTOMATIC
  integer :: istate

  if (velflag.ne.2) then
     OFLWR "what? velflag check velocityop",velflag; CFLST
  endif

  gaugetransform(:,:)= (0d0,-1d0) * avector * couplingmatlen(:,:)

  call expmat(gaugetransform,numstates)

  h0mat=0d0
  do istate=1,numstates
     h0mat(istate,istate)=myStateEnergies(istate)
  enddo

  velocityop(:,:) = MATMUL(gaugetransform(:,:),MATMUL(h0mat(:,:),CONJG(gaugetransform(:,:)))) &
       - h0mat(:,:)

end subroutine get_velocityop





