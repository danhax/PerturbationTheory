
#include "Definitions.INC"


program pert
  use parameters
  implicit none
  integer :: i

  pi=4d0*atan(1d0)
  do i=1,SLN
     nullbuff(i:i)=" "
  enddo

  call getinpfile()
  call getparams()
  call perttheory()

end program pert


subroutine perttheory()
  use parameters
  use pulsesubmod
  implicit none
  real*8 :: t1,t2
  integer :: itime,numiters
  complex*16, allocatable :: zerovec(:,:), onevec(:,:), twovec(:,:), threevec(:,:), sumvec(:,:)

  allocate(zerovec(numstates,0:numsteps), onevec(numstates,0:numsteps), &
       twovec(numstates,0:numsteps), threevec(numstates,0:numsteps), sumvec(numstates,0:numsteps))

!! zeroth-order wave function
  zerovec(:,:)=0d0
  do itime=0,numsteps
     zerovec(1,itime)=exp((0d0,-1d0)*stateEnergies(1)*itime*par_timestep)
  enddo

!!!!!    first-order wave function    !!!!!

  onevec(:,:)=0d0

  do itime=1,numsteps  !! propagation from itime-1 to itime
     t1=(itime-1)*par_timestep
     t2=itime*par_timestep

     call expoprop(t1,zerovec(:,itime-1),t2,zerovec(:,itime),onevec(:,itime-1),onevec(:,itime),numiters)
  enddo

  open(1111, file="onevec.dat",status="unknown")
  do itime=0,numsteps
     write(1111,'(1000F18.10)') itime*par_timestep, vectdpot(itime*par_timestep,velflag), onevec(:,itime)
  enddo
  close(1111)
  open(1111, file="onevec_abs.dat",status="unknown")
  do itime=0,numsteps
     write(1111,'(1000F18.10)') itime*par_timestep, vectdpot(itime*par_timestep,velflag), abs(onevec(:,itime)**2)
  enddo
  close(1111)

!!!!!!!   second-order wave function   !!!!!!

  twovec(:,:)=0d0

  do itime=1,numsteps  !! propagation from itime-1 to itime
     t1=(itime-1)*par_timestep
     t2=itime*par_timestep

     call expoprop(t1,onevec(:,itime-1),t2,onevec(:,itime),twovec(:,itime-1),twovec(:,itime),numiters)
  enddo

  open(1111, file="twovec.dat",status="unknown")
  do itime=0,numsteps
     write(1111,'(1000F18.10)') itime*par_timestep, vectdpot(itime*par_timestep,velflag), twovec(:,itime)
  enddo
  close(1111)
  open(1111, file="twovec_abs.dat",status="unknown")
  do itime=0,numsteps
     write(1111,'(1000F18.10)') itime*par_timestep, vectdpot(itime*par_timestep,velflag), abs(twovec(:,itime)**2)
  enddo
  close(1111)


!!!!!!!   third-order wave function   !!!!!!

  threevec(:,:)=0d0

  do itime=1,numsteps  !! propagation from itime-1 to itime
     t1=(itime-1)*par_timestep
     t2=itime*par_timestep

     call expoprop(t1,twovec(:,itime-1),t2,twovec(:,itime),threevec(:,itime-1),threevec(:,itime),numiters)
  enddo

  open(1111, file="threevec.dat",status="unknown")
  do itime=0,numsteps
     write(1111,'(1000F18.10)') itime*par_timestep, vectdpot(itime*par_timestep,velflag), threevec(:,itime)
  enddo
  close(1111)
  open(1111, file="threevec_abs.dat",status="unknown")
  do itime=0,numsteps
     write(1111,'(1000F18.10)') itime*par_timestep, vectdpot(itime*par_timestep,velflag), abs(threevec(:,itime)**2)
  enddo
  close(1111)


!! SUM THEM

!!  sumvec(:,:)=zerovec(:,:) + onevec(:,:) + twovec(:,:)

  sumvec(:,:)=zerovec(:,:) + onevec(:,:) + twovec(:,:) + threevec(:,:)

  open(1111, file="sumvec.dat",status="unknown")
  do itime=0,numsteps
     write(1111,'(1000F18.10)') itime*par_timestep, vectdpot(itime*par_timestep,velflag), sumvec(:,itime)
  enddo
  close(1111)
  open(1111, file="sumvec_abs.dat",status="unknown")
  do itime=0,numsteps
     write(1111,'(1000F18.10)') itime*par_timestep, vectdpot(itime*par_timestep,velflag), abs(sumvec(:,itime)**2), &
          SUM(abs(sumvec(:,itime)**2))
  enddo
  close(1111)

  deallocate(zerovec, onevec, twovec, threevec, sumvec)

  OFLWR "done perttheory"; CFLST

end subroutine perttheory



module jacopmod

contains
  subroutine jacoperate(invec,outvec)
    use parameters
    implicit none
    complex*16,intent(in) :: invec(numstates)
    complex*16,intent(out) :: outvec(numstates)
    integer :: istate

    do istate=1,numstates
       outvec(istate)=invec(istate)*(0d0,-1d0)*stateEnergies(istate)
    enddo
  end subroutine jacoperate

  subroutine potoperate(midtime,invec,outvec)
    use parameters
    use pulsesubmod
    implicit none
    real*8, intent(in) :: midtime
    complex*16, intent(in) :: invec(numstates)
    complex*16, intent(out) :: outvec(numstates)
    complex*16 :: mypot

    mypot = vectdpot(midtime, velflag)
    outvec(:) = mypot * MATMUL(couplingmat(1:numstates,1:numstates), invec(:))
  end subroutine potoperate

end module jacopmod



subroutine expoprop(time1, drivingvec1, time2, drivingvec2, invec, outvec, numiters)
  use parameters
  use jacopmod
  implicit none
  real*8,intent(in) :: time1,time2
  complex*16,intent(in) :: invec(numstates), drivingvec1(numstates), drivingvec2(numstates)
  complex*16,intent(out) :: outvec(numstates)
  integer,intent(out) :: numiters
  complex*16 :: avec1(numstates), avec2(numstates)
  real*8 :: tdiff, midtime
  integer, save :: icalled=0
  integer :: istate

  numiters=0
  icalled=icalled+1

  avec1=0; avec2=0
  tdiff=(time2-time1)
  midtime=(time2+time1)/2d0

  call potoperate(midtime, drivingvec1(:), avec1(:))
  call potoperate(midtime, drivingvec2(:), avec2(:))
  avec1(:)=avec1(:)*(0d0,-1d0)
  avec2(:)=avec2(:)*(0d0,-1d0)

  do istate=1,numstates

!! zeroth order
     outvec(istate)=exp((0d0,-1d0)*stateEnergies(istate)*par_timestep)*invec(istate) + &
          par_timestep*myphi1((0d0,-1d0)*stateEnergies(istate)*par_timestep)*(avec1(istate)+avec2(istate))/2d0

!! first order... think there is a mistake, not working better
!     outvec(istate)=exp((0d0,-1d0)*stateEnergies(istate)*par_timestep)*invec(istate) + &
!          par_timestep*myphi1((0d0,-1d0)*stateEnergies(istate)*par_timestep)*avec1(istate) + &
!          par_timestep*myphi2((0d0,-1d0)*stateEnergies(istate)*par_timestep)*(avec2(istate)-avec1(istate))

  enddo

contains
  function myphi1(zval)
    complex*16 :: myphi1
    complex*16, intent(in) :: zval
    real*8, parameter :: minval=1d-6
    if (abs(zval).lt.minval) then
       myphi1 = 1 + zval/2 + zval**2/6
    else
       myphi1 = (exp(zval)-1)/zval
    endif
  end function myphi1

  function myphi2(zval)
    complex*16 :: myphi2
    complex*16, intent(in) :: zval
    real*8, parameter :: minval=1d-6
    if (abs(zval).lt.minval) then
       myphi2 = 0.5d0 + zval/6 + zval**2/24
    else
       myphi2 = (myphi1(zval)-1)/zval
    endif
  end function myphi2

end subroutine expoprop


subroutine mpistop()
  use fileptrmod
  implicit none
  call waitawhile()
  OFLWR "PERTTHEORY STOP!";CFL
  stop
end subroutine mpistop

subroutine openfile()
end subroutine openfile

subroutine closefile()
end subroutine closefile

