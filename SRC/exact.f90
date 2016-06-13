
#include "Definitions.INC"


subroutine schrodinger()
  use parameters
  use pulsesubmod
  implicit none
  real*8 :: t1,t2
  integer :: itime
  complex*16, allocatable :: propvec(:,:)

  allocate(propvec(numstates,0:numsteps))

  propvec(:,:)=0d0
  propvec(1,0) = 1d0

  do itime=1,numsteps  !! propagation from itime-1 to itime
     t1=(itime-1)*par_timestep
     t2=itime*par_timestep
     call exactprop(t1,t2,propvec(:,itime-1),propvec(:,itime))
  enddo

  call save1(numstates,0,numsteps,par_timestep,propvec(:,:),"propvec",0)
  call save1(numstates,0,numsteps,par_timestep,propvec(:,:),"propvec",1)

  deallocate(propvec)

  OFLWR "done schrodinger"; CFLST

end subroutine schrodinger


module exactmod
  implicit none
  complex*16, allocatable :: hammatexp(:,:),h0diagexp(:)
end module exactmod

subroutine exactalloc()
  use parameters
  use exactmod
  implicit none
  allocate(hammatexp(numstates,numstates),h0diagexp(numstates))
  hammatexp=0d0
  h0diagexp(:) = exp((0d0,-1d0) * par_timestep * myStateEnergies(1:numstates))
end subroutine exactalloc


subroutine exactprop(time1,time2,invec,outvec)
  use parameters
  use pulsesubmod
  use exactmod
  use expsubmod
  implicit none
  real*8,intent(in) :: time1,time2
  complex*16, intent(in) :: invec(numstates)
  complex*16, intent(out) :: outvec(numstates)
  real*8 :: mypot,midtime,tdiff
  integer :: istate

  tdiff=(time2-time1)
  midtime=(time2+time1)/2d0
  mypot=vectdpot(midtime, velflag)

  if (mypot.eq.0d0) then
     outvec(:) = invec(:) * h0diagexp(:)
  else
     hammatexp(:,:) = mypot * (0d0,-1d0) * par_timestep * mycouplingmat(1:numstates,1:numstates)
     do istate=1,numstates
        hammatexp(istate,istate) = hammatexp(istate,istate) + &
             (0d0,-1d0) * par_timestep * myStateEnergies(istate)
     enddo
     call expmat(hammatexp,numstates)
     outvec(:) = MATMUL(hammatexp(:,:),invec(:))
  endif

end subroutine exactprop


