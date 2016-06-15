
#include "Definitions.INC"


subroutine perttheory()
  use parameters
  use pulsesubmod
  implicit none
  real*8 :: t1,t2
  integer :: itime
  complex*16, allocatable :: zerovec(:,:), onevec(:,:), twovec(:,:), threevec(:,:), &
       fourvec(:,:), sumvec(:,:)

  allocate(zerovec(numstates,0:numsteps), onevec(numstates,0:numsteps), &
       twovec(numstates,0:numsteps), threevec(numstates,0:numsteps), &
       fourvec(numstates,0:numsteps), sumvec(numstates,0:numsteps))

!! zeroth-order wave function
  zerovec(:,:)=0d0
  do itime=0,numsteps
     zerovec(1,itime)=exp((0d0,-1d0)*myStateEnergies(1)*itime*par_timestep)
  enddo

!!!!!    first-order wave function    !!!!!

  onevec(:,:)=0d0

  if (maxorder.ge.1) then

  do itime=1,numsteps  !! propagation from itime-1 to itime
     t1=(itime-1)*par_timestep
     t2=itime*par_timestep
     call expoprop(t1,zerovec(:,itime-1),t2,zerovec(:,itime),onevec(:,itime-1),onevec(:,itime))
  enddo

  call save1(0,numsteps,par_timestep,onevec(:,:),"onevec",0)
  call save1(0,numsteps,par_timestep,onevec(:,:),"onevec",1)

  endif

!!!!!!!   second-order wave function   !!!!!!

  twovec(:,:)=0d0

  if (maxorder.ge.2) then

  do itime=1,numsteps  !! propagation from itime-1 to itime
     t1=(itime-1)*par_timestep
     t2=itime*par_timestep
     call expoprop(t1,onevec(:,itime-1),t2,onevec(:,itime),twovec(:,itime-1),twovec(:,itime))
  enddo

  call save1(0,numsteps,par_timestep,twovec(:,:),"twovec",0)
  call save1(0,numsteps,par_timestep,twovec(:,:),"twovec",1)

  endif

!!!!!!!   third-order wave function   !!!!!!

  threevec(:,:)=0d0

  if (maxorder.ge.3) then

  do itime=1,numsteps  !! propagation from itime-1 to itime
     t1=(itime-1)*par_timestep
     t2=itime*par_timestep

     call expoprop(t1,twovec(:,itime-1),t2,twovec(:,itime),threevec(:,itime-1),threevec(:,itime))
  enddo

  call save1(0,numsteps,par_timestep,threevec(:,:),"threevec",0)
  call save1(0,numsteps,par_timestep,threevec(:,:),"threevec",1)

  endif

!!!!!!!   fourth-order wave function   !!!!!!

  fourvec(:,:)=0d0

  if (maxorder.ge.4) then

  do itime=1,numsteps  !! propagation from itime-1 to itime
     t1=(itime-1)*par_timestep
     t2=itime*par_timestep

     call expoprop(t1,threevec(:,itime-1),t2,threevec(:,itime),fourvec(:,itime-1),fourvec(:,itime))
  enddo

!  call save1(0,numsteps,par_timestep,fourvec(:,:),"fourvec",0)
!  call save1(0,numsteps,par_timestep,fourvec(:,:),"fourvec",1)

  endif

  if (maxorder.gt.4) then
     OFLWR "maxorder.gt.4 not supported",maxorder; CFLST
  endif

!! SUM THEM

  sumvec(:,:)=zerovec(:,:) + onevec(:,:) + twovec(:,:) + threevec(:,:) + fourvec(:,:)

  call save1(0,numsteps,par_timestep,sumvec(:,:),"sumvec",0)
  call save1(0,numsteps,par_timestep,sumvec(:,:),"sumvec",1)

  deallocate(zerovec, onevec, twovec, threevec, fourvec, sumvec)

  OFLWR "done perttheory"; CFLST

end subroutine perttheory



module jacopmod
  implicit none

contains
  subroutine jacoperate(invec,outvec)
    use parameters
    implicit none
    complex*16,intent(in) :: invec(numstates)
    complex*16,intent(out) :: outvec(numstates)
    integer :: istate

    do istate=1,numstates
       outvec(istate)=invec(istate)*(0d0,-1d0)*myStateEnergies(istate)
    enddo
  end subroutine jacoperate

  subroutine potoperate(midtime,invec,outvec)
    use parameters
    use pulsesubmod
    implicit none
    real*8, intent(in) :: midtime
    complex*16, intent(in) :: invec(numstates)
    complex*16, intent(out) :: outvec(numstates)
    complex*16 :: velocityop(numstates,numstates) !! AUTOMATIC
    real*8 :: mypot

    velocityop(:,:)=0d0

    mypot = vectdpot(midtime, velflag)

    if (velflag.eq.2) then
       outvec(:)=0d0
    else
       outvec(:) = mypot * MATMUL(mycouplingmat(1:numstates,1:numstates), invec(:))
    endif

    if (velflag.eq.1) then
!! factor of 1/4 in asquaredop
       outvec(:) = outvec(:) + mypot**2 * MATMUL(asquaredop(:,:),invec(:))
    elseif (velflag.eq.2) then
       call get_velocityop(mypot,velocityop(:,:))
       outvec(:) = outvec(:) + MATMUL(velocityop(:,:),invec(:))
    endif

  end subroutine potoperate

end module jacopmod



subroutine expoprop(time1, drivingvec1, time2, drivingvec2, invec, outvec)
  use parameters
  use jacopmod
  implicit none
  real*8,intent(in) :: time1,time2
  complex*16,intent(in) :: invec(numstates), drivingvec1(numstates), drivingvec2(numstates)
  complex*16,intent(out) :: outvec(numstates)
  complex*16 :: avec1(numstates), avec2(numstates)
  real*8 :: tdiff, midtime
  integer, save :: icalled=0
  integer :: istate

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
     outvec(istate)=exp((0d0,-1d0)*myStateEnergies(istate)*par_timestep)*invec(istate) + &
          par_timestep*myphi1((0d0,-1d0)*myStateEnergies(istate)*par_timestep)*(avec1(istate)+avec2(istate))/2d0

!! first order... think there is a mistake, not working better
!     outvec(istate)=exp((0d0,-1d0)*myStateEnergies(istate)*par_timestep)*invec(istate) + &
!          par_timestep*myphi1((0d0,-1d0)*myStateEnergies(istate)*par_timestep)*avec1(istate) + &
!          par_timestep*myphi2((0d0,-1d0)*myStateEnergies(istate)*par_timestep)*(avec2(istate)-avec1(istate))

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

