
!!
!!  UTILITIES.
!!  

#include "Definitions.INC"

#ifdef BOOGABOOGAAAAH

subroutine complexdiff(size,in,out,docirc)
  implicit none
  integer, intent(in) :: size
  logical, intent(in) :: docirc
  complex*16, intent(in) :: in(size)
  complex*16, intent(out) :: out(size)
  integer :: i,jj

  if (docirc) then
!! guarantees F.T. at zero is zero, right?

     do i=1,size
        out(i)= &
             1d0/280d0 * in(cindex(i-4)) &
             - 4d0/105d0 * in(cindex(i-3)) &
             + 1d0/5d0 * in(cindex(i-2)) &
             - 4d0/5d0 * in(cindex(i-1)) &
             + 4d0/5d0 * in(cindex(i+1)) &
             - 1d0/5d0 * in(cindex(i+2)) &
             + 4d0/105d0 * in(cindex(i+3)) &
             - 1d0/280d0 * in(cindex(i+4))
     enddo

  else

     out(:)=0d0
     do i=2,size-1
        jj=min(min(i-1,size-i),4)
        select case(jj)
        case(1)
           out(i)= &
                - 1d0/2d0 * in(i-1) &
                + 1d0/2d0 * in(i+1)
        case(2)
           out(i)= &
                1d0/12d0 * in(i-2) &
                - 2d0/3d0 * in(i-1) &
                + 2d0/3d0 * in(i+1) &
                - 1d0/12d0 * in(i+2)
        case(3)
           out(i)= &
                - 1d0/60d0 * in(i-3) &
                + 3d0/20d0 * in(i-2) &
                - 3d0/4d0 * in(i-1) &
                + 3d0/4d0 * in(i+1) &
                - 3d0/20d0 * in(i+2) &
                + 1d0/60d0 * in(i+3)
        case(4)
           out(i)= &
                1d0/280d0 * in(i-4) &
                - 4d0/105d0 * in(i-3) &
                + 1d0/5d0 * in(i-2) &
                - 4d0/5d0 * in(i-1) &
                + 4d0/5d0 * in(i+1) &
                - 1d0/5d0 * in(i+2) &
                + 4d0/105d0 * in(i+3) &
                - 1d0/280d0 * in(i+4)
        end select
     end do

  endif  !! docirc

contains
  function cindex(inindex)
    integer :: cindex,inindex
    cindex=mod(2*size+inindex-1,size)+1
  end function cindex

end subroutine complexdiff


subroutine zfftf_wrap_diff(size,inout,diffdflag)
  implicit none
  integer, intent(in) :: size,diffdflag
  complex*16, intent(inout) :: inout(size)
  complex*16,allocatable :: work(:)
  integer :: i

  if (diffdflag.eq.0) then
     call zfftf_wrap(size,inout)
  else
     allocate(work(size)); work=0

!! guarantees F.T. at zero is zero, right?
#define DOCIRC .true.

     call complexdiff(size,inout,work,DOCIRC)

     call zfftf_wrap(size,work)

     do i=1,size
        inout(i)=work(i)*facfunct(i-1,size-1,1)
     enddo

     deallocate(work)
  endif

contains
  function facfunct(myindex,numdata,diffdflag)
    use fileptrmod
    implicit none
    integer, intent(in) :: myindex,diffdflag,numdata
    complex*16 :: facfunct,ccsum
    real*8, parameter :: twopi = 6.28318530717958647688d0
    if (myindex.lt.0.or.myindex.gt.numdata) then
       OFLWR "FACFUNCT ERR", myindex,0,numdata; CFLST
    endif
    ccsum=1d0
    if (diffdflag.ne.0) then
       if (myindex.ne.0) then
          ccsum= 1d0 / ((0d0,1d0)*myindex) / twopi * (numdata+1)
       else
          ccsum=0d0
       endif
    endif
    facfunct=ccsum
  end function facfunct

end subroutine zfftf_wrap_diff


subroutine zfftf_wrap(size,inout)
  implicit none
  integer, intent(in) :: size
  complex*16, intent(inout) :: inout(size)
  complex*16,allocatable :: wsave(:)
  allocate(wsave(4*size+15))
  wsave(:)=0d0
  call zffti(size,wsave)
  call zfftf(size,inout,wsave)
  deallocate(wsave)
end subroutine zfftf_wrap


subroutine zfftb_wrap(size,inout)
  implicit none
  integer, intent(in) :: size
  complex*16, intent(inout) :: inout(size)
  complex*16,allocatable :: wsave(:)
  allocate(wsave(4*size+15))
  wsave(:)=0d0
  call zffti(size,wsave)
  call zfftb(size,inout,wsave)
  deallocate(wsave)
end subroutine zfftb_wrap

#endif

subroutine waitawhile()
  implicit none
  integer :: i,j
  character (len=10) :: mytext
  j=1
  do i=1,10000000
     j=mod(j*i,1777)
  enddo
  write(mytext,'(I10)') j
  call system("echo "//mytext//" >> /dev/null")
end subroutine waitawhile


subroutine checkiostat(iniostat,intext)
  use fileptrmod
  implicit none
  integer,intent(in) :: iniostat
  character*(*),intent(in) :: intext
  if (iniostat /=0 ) then
     print *, "ABORT: I/O error ", iniostat,intext
     OFLWR "ABORT: I/O error ", iniostat,intext; CFL
     call waitawhile()
     stop
     stop   !!   STOP.   !!
     stop
  endif
end subroutine checkiostat

!! v1.27 getlen now reports length of string not length of string plus one

function getlen(buffer)
  implicit none
  character buffer*(*)
  integer :: j, getlen, mylen
  mylen=LEN(buffer)
  j=1
  do while (j.le.mylen)
     if (buffer(j:j) .eq. " ") then
        getlen=j-1
        return
     else
        j=j+1
     endif
  enddo
  getlen=mylen
end function getlen


function getlen2(buffer)
  implicit none
  character buffer*(*)
  integer :: j, getlen2, nn
  nn=LEN(buffer)-4
  j=1
  do while ((j.lt.nn).and..not.(buffer(j:j+3) .eq. "    "))
     j=j+1
  enddo
  getlen2=j-1
end function getlen2


function floatfac(in)
  implicit none
  integer,intent(in) :: in
  integer ::  i
  real*8 :: floatfac, sum
  sum=1.d0
  do i=1,in
     sum=sum*i
  enddo
  floatfac=sum
end function floatfac


function myisnan(input)
  implicit none
  real*8,intent(in) :: input
  logical :: myisnan
  if ((input+1.0d0.eq.input)) then
     myisnan=.true.
  else
     myisnan=.false.
  endif
end function myisnan




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      MODULE INVSUBMOD    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module invsubmod
contains

!!$  if relflag then
!!$  tol is the maximum ratio of smallest to biggest eigenvalue
!!$    i.e. tol is relative not absolute
!!$  it .not.relflag then tol is absolute

subroutine invmatsmooth(A,N,LDA,tol,relflag)  !! inverse of ANY matrix.
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=(A_in)**-1
  implicit none
  logical,intent(in) :: relflag
  integer,intent(in) :: N,LDA
  real*8,intent(in) :: tol
  real*8 :: mytol
  complex*16,intent(inout) :: A(LDA,N)
  integer :: lwork,i,j,k
  real*8,allocatable :: SV(:)
  complex*16,allocatable :: SAVEA(:,:),U(:,:),VT(:,:),work(:)
  complex*16,allocatable :: zwork(:)

  allocate(zwork(5*N))
  allocate(SV(N),SAVEA(LDA,N),U(N,N),VT(N,N),work(5*N))

  lwork=5*N

  SAVEA=A

  call zgesvd('A','A',N,N,A,LDA,SV,U,N,VT,N,zwork,lwork,work,i)

  if (i.ne.0) then
     print *, "ERR SVD",i;
     do i=1,N
        print *, SAVEA(:,i)
     enddo
     stop
  endif

  if (relflag) then 
     mytol=tol*SV(1)
  else
     mytol=tol
  endif

!! apply the function
  do k=1,N
     if (SV(k).ne.0d0) then
        if(SV(k).lt.mytol) then    !! it is positive ,whatever

!!$ NAAAH    SV(k)= 1d0 / mytol * (3 * (SV(k) / mytol / ) - 2 *( SV(k) / mytol )**2)

!!$    SVD so SV is real, keeping old regularization for now v1.19

           SV(k)= 1d0 / mytol

        else
           SV(k) = 1d0 / SV(k)
        endif
     endif
  enddo

!! rebuild the matrix
  do j=1,N
    do i=1,N
      A(i,j) = 0d0
      do k=1,N
         A(i,j) = A(i,j) + CONJG(VT(k,i)) * SV(k) * CONJG(U(j,k))
       enddo
    enddo
  enddo

  deallocate(zwork)
  deallocate(SV,SAVEA,U,VT,work)

end subroutine invmatsmooth



!!$  if relflag then
!!$  tol is the maximum ratio of smallest to biggest eigenvalue
!!$    i.e. tol is relative not absolute
!!$  it .not.relflag then tol is absolute


subroutine realinvmatsmooth(A,N,tol,relflag)  !! inverse of ANY matrix.
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=(A_in)**-1
  implicit none
  logical,intent(in) :: relflag
  integer,intent(in) :: N
  real*8,intent(in) :: tol 
  real*8 :: mytol
  real*8,intent(inout) :: A(N,N)
  real*8,allocatable :: U(:,:),VT(:,:),work(:),SV(:)
  integer :: lwork,i,j,k

  allocate( U(N,N),VT(N,N),work(5*N),SV(N) )

  lwork=5*N

!! do the svd

  call dgesvd('A','A',N,N,A,N,SV,U,N,VT,N,work,lwork,i)
  if (i.ne.0) then
     print *, "ERR SVD";     stop
  endif

  if (relflag) then 
     mytol=tol*SV(1)
  else
     mytol=tol
  endif

!! apply the function
  do k=1,N
     if (SV(k).ne.0d0) then
        if(abs(SV(k)).lt.mytol) then

!!$ NAAAH         SV(k)= 1d0 / mytol * (3 * (SV(k) / mytol) - 2 *( SV(k) / mytol )**2)

!!$    SVD so SV is real, keeping old regularization for now v1.19

           SV(k) = 1d0 / mytol 
        else
           SV(k) = 1d0 / SV(k)
        endif
     endif
  enddo
!! rebuild the matrix
  do j=1,N
    do i=1,N
      A(i,j) = 0d0
      do k=1,N
         A(i,j) = A(i,j) + VT(k,i) * SV(k) * U(j,k)
       enddo
    enddo
  enddo

  deallocate( U,VT,work,SV)

end subroutine realinvmatsmooth

end module invsubmod




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      MODULE EXPSUBMOD    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module expsubmod
contains

subroutine expmat(A,N)
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=exp(A_in)
  use tol_parameters
  use invsubmod
  implicit none
  integer,intent(in) :: N
  complex*16,intent(inout) :: A(N,N)
  integer :: lwork,i,k,j
  complex*16,allocatable :: eig(:), CVL(:,:),CVR(:,:),CA(:,:)
  complex*16,allocatable :: VL(:,:),VR(:,:),work(:),rwork(:)

  allocate( eig(N), CVL(N,N),CVR(N,N),CA(N,N), VL(N,N),VR(N,N),work(8*N),rwork(4*N) )
  lwork=8*N

  j=0
  call zgeev('V','V',N,A,N,eig,VL,N,VR,N,work,lwork,rwork,i)

  VL(:,:)=TRANSPOSE(VR(:,:))
  call invmatsmooth(VL,N,N,invtol,.true.)

!! apply the function
  do k=1,N
     eig(k) = exp( eig(k) )
     CVR(:,k)=VR(:,k)*eig(k)
  enddo
  CVL(:,:)=VL(:,:)

!! rebuild the matrix

  call ZGEMM('N','T',N,N,N,(1d0,0d0),CVR,N,CVL,N,(0d0,0d0),CA,N)
  A(:,:)=CA(:,:)  !! OK IMP CONV

  deallocate( eig,cvl,cvr,ca,vl,vr,work,rwork )

end subroutine expmat

end module expsubmod


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      MODULE MATSUBMOD    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! orthog routines:
!!  FLAG 1:  if on input A(i,j) = <phi_i|phi_j>
!!              then on output varphi_j = sum_i A(i,j) |phi_i>  are orthonormal.
!!  FLAG 2:  if on input A(i,j) = <phi_i|phi_j>
!!              on output varphi_j = sum_i A(i,j) |phi_i> are biorthonormal to phi:
!!            : <varphi_i|phi_j>=delta_ij


module matsubmod
contains

subroutine allpurposemat(A,N,flag)
  use tol_parameters
  use fileptrmod
  use invsubmod
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=sqrt^-1(A_in) if flag=1 or sqrt(A_in) if flag=2
  implicit none
  integer,intent(in) :: N,flag
  complex*16,intent(inout) :: A(N,N)
  complex*16,allocatable :: eig(:), CVL(:,:),CVR(:,:),CA(:,:)
  complex*16,allocatable :: VL(:,:),VR(:,:),work(:),rwork(:)
  integer :: lwork,i,k,j

  allocate(eig(N), CVL(N,N),CVR(N,N),CA(N,N), VL(N,N),VR(N,N),work(8*N),rwork(4*N))
  lwork=8*N

  j=0
  call zgeev('V','V',N,A,N,eig,VL,N,VR,N,work,lwork,rwork,i)

  VL(:,:)=TRANSPOSE(VR(:,:))
  call invmatsmooth(VL,N,N,invtol,.true.)

!! apply the function
  do k=1,N
     select case(flag)
     case(1)
        eig(k) = 1d0 / sqrt( eig(k) )
     case(2)
        eig(k) = sqrt( eig(k) )
     case default
        OFLWR "OOGA ALLPURPOSEMAT ",flag; CFLST
     end select
     CVR(:,k)=VR(:,k)*eig(k)
  enddo
  CVL(:,:)=VL(:,:)

!! rebuild the matrix

  call ZGEMM('N','T',N,N,N,(1d0,0d0),CVR,N,CVL,N,(0d0,0d0),CA,N)
  A(:,:)=CA(:,:)  !! OK IMP CONV

  deallocate(eig, CVL,CVR,CA, VL,VR,work,rwork)

end subroutine allpurposemat

end module matsubmod



subroutine assigncomplex(realmat,complexf)
  implicit none
  complex*16,intent(in) :: complexf
  real*8,intent(out) :: realmat(2,2)
  realmat(1,1)=real(complexf,8);  realmat(2,2)=real(complexf,8)
  realmat(2,1)=imag(complexf);  realmat(1,2)=(-1)*imag(complexf)
end subroutine assigncomplex

subroutine assigncomplexmat(realmat,complexf,m,n)
  implicit none
  integer,intent(in) :: n,m
  complex*16,intent(in) :: complexf(m,n)
  real*8,intent(out) :: realmat(2,m,2,n)
  realmat(1,:,1,:)=real(complexf(:,:),8);  realmat(2,:,2,:)=real(complexf(:,:),8)
  realmat(2,:,1,:)=imag(complexf(:,:));  realmat(1,:,2,:)=(-1)*imag(complexf(:,:))
end subroutine assigncomplexmat

subroutine assigncomplexvec(realmat,complexf,m)
  implicit none
  integer,intent(in) :: m
  complex*16,intent(in) :: complexf(m)
  real*8,intent(out) :: realmat(2,m)
  realmat(1,:)=real(complexf(:),8);  realmat(2,:)=imag(complexf(:))
end subroutine assigncomplexvec

subroutine assignrealvec(complexf,realmat,m)
  implicit none
  integer,intent(in) :: m
  complex*16,intent(out) :: complexf(m)
  real*8,intent(in) :: realmat(2,m)
  complexf(:)=realmat(1,:)+realmat(2,:)*(0d0,1d0)
end subroutine assignrealvec


