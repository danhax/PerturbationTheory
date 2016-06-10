
#include "Definitions.INC"

module pulsesubmod
contains

function vectdpot(myintime,invelflag)
  use pulse_parameters
  implicit none
  real*8,intent(in) :: myintime
  integer, intent(in) :: invelflag
  real*8 :: vectdpot, tdpotout

  call vectdpot0(myintime,invelflag,tdpotout,1,numpulses)
  vectdpot=tdpotout

end function vectdpot


subroutine vectdpot0(myintime,invelflag,tdpotout,ilow,ihigh)
  use pulse_parameters
  use fileptrmod
  implicit none
  integer, intent(in) :: invelflag,ilow,ihigh
  real*8,intent(in) :: myintime
  real*8,intent(out) :: tdpotout

  if (invelflag.eq.0) then
     tdpotout=tdpotlen0(myintime,ilow,ihigh)
  else
     tdpotout=tdpotvel0(myintime,ilow,ihigh)
  endif

end subroutine vectdpot0


function tdpotlen0(myintime, ilow,ihigh)
  use pulse_parameters
  use fileptrmod
  implicit none
  integer,intent(in) :: ilow,ihigh
  real*8,intent(in) :: myintime
  integer :: ipulse
  real*8 :: tdpotlen0

  tdpotlen0=0.d0

!!  do ipulse=1,numpulses

  do ipulse=ilow,ihigh
     select case (pulsetype(ipulse))
     case (1)
        tdpotlen0=tdpotlen0+simplepulselen(myintime,ipulse)
     case (2)
        tdpotlen0=tdpotlen0+pulselen(myintime,ipulse)
     case (3)
        tdpotlen0=tdpotlen0+longpulselen(myintime,ipulse)
     case (4)
        tdpotlen0=tdpotlen0+cwpulselen(myintime,ipulse)
     case default
        OFLWR "Pulse type not supported: ", pulsetype(ipulse); CFLST
     end select
  enddo
  
end function tdpotlen0


function tdpotvel0(myintime,ilow,ihigh)
  use pulse_parameters
  use fileptrmod
  implicit none
  integer,intent(in) :: ilow,ihigh
  real*8,intent(in) :: myintime
  integer :: ipulse
  real*8 :: tdpotvel0

  tdpotvel0=0.d0

!!  do ipulse=1,numpulses

  do ipulse=ilow,ihigh
     select case (pulsetype(ipulse))
     case (1)
        tdpotvel0=tdpotvel0+simplepulsevel(myintime, ipulse)
     case (2)
        tdpotvel0=tdpotvel0+pulsevel(myintime, ipulse)
     case (3)
        tdpotvel0=tdpotvel0+longpulsevel(myintime, ipulse)
     case (4)
        tdpotvel0=tdpotvel0+cwpulsevel(myintime, ipulse)
     case default
        OFLWR "Pulse type not supported: ", pulsetype(ipulse); CFLST
     end select
  enddo

end function tdpotvel0


function simplepulselen(myintime, ipulse)
  use pulse_parameters
  use constant_parameters
  implicit none
  integer,intent(in) :: ipulse
  real*8,intent(in) :: myintime
  real*8 :: time
  real*8 :: simplepulselen

  simplepulselen=0.d0

  if (myintime.ge.pulsestart(ipulse)) then
     time=myintime-pulsestart(ipulse)
     if (time.le.pi/omega(ipulse)) then

        simplepulselen = pulsestrength(ipulse) * omega(ipulse) * &
             ( sin(time*omega(ipulse))*cos(time*omega(ipulse) + phaseshift(ipulse)) &
             + sin(time*omega(ipulse) + phaseshift(ipulse))*cos(time*omega(ipulse)) )

     endif
  endif
end function simplepulselen


function simplepulsevel(myintime, ipulse)
  use pulse_parameters
  use constant_parameters
  implicit none
  integer,intent(in) :: ipulse
  real*8,intent(in) :: myintime
  real*8 :: time
  real*8 :: simplepulsevel

  simplepulsevel=0.d0

  if (myintime.ge.pulsestart(ipulse)) then
     time=myintime-pulsestart(ipulse)
     if (time.le.pi/omega(ipulse)) then

        simplepulsevel = pulsestrength(ipulse) * &
             sin(time*omega(ipulse)) * sin(time*omega(ipulse) + phaseshift(ipulse) )

     endif
  endif

end function simplepulsevel


function cwpulselen(myintime, ipulse)
  use pulse_parameters
  use constant_parameters
  implicit none
  integer,intent(in) :: ipulse
  real*8,intent(in) :: myintime
  real*8 :: cwpulselen

  cwpulselen=0d0
  if (myintime.lt.pi/omega(ipulse)) then
     cwpulselen = pulsestrength(ipulse) * omega2(ipulse) * &
          cos(myintime*omega2(ipulse)+phaseshift(ipulse))
  endif

end function cwpulselen


function cwpulsevel(myintime, ipulse)
  use pulse_parameters
  use constant_parameters
  implicit none
  integer,intent(in) :: ipulse
  real*8,intent(in) :: myintime
  real*8 :: cwpulsevel

  cwpulsevel=0.d0
  if (myintime.lt.pi/omega(ipulse)) then
     cwpulsevel = pulsestrength(ipulse) * sin(myintime*omega2(ipulse) + phaseshift(ipulse))
  endif

end function cwpulsevel


function pulselen(myintime, ipulse)
  use pulse_parameters
  use constant_parameters
  use fileptrmod
  implicit none
  integer,intent(in) :: ipulse
  real*8,intent(in) :: myintime
  real*8 :: time
  real*8 :: pulselen

  pulselen=0.d0

  if (chirp(ipulse).ne.0d0.or.ramp(ipulse).ne.0d0) then
     OFLWR "Chirp / ramp not supported for length", chirp(ipulse), ipulse; CFLST
  endif
  if (myintime.ge.pulsestart(ipulse)) then
     time=myintime-pulsestart(ipulse)
     if (time.le.pi/omega(ipulse)) then
        pulselen = pulsestrength(ipulse) * ( &
             2*omega(ipulse)*sin(time*omega(ipulse))*cos(time*omega(ipulse)) * &
             sin((time-pi/omega(ipulse)/2)*omega2(ipulse) + phaseshift(ipulse)) &
             + sin(time*omega(ipulse))**2 * omega2(ipulse) * &
             cos((time-pi/omega(ipulse)/2)*omega2(ipulse) + phaseshift(ipulse)) )
     endif
  endif

end function pulselen


function pulsevel(myintime, ipulse)
  use pulse_parameters
  use constant_parameters
  implicit none
  integer,intent(in) :: ipulse
  real*8,intent(in) :: myintime
  real*8 :: time,thisomega2
  real*8 :: pulsevel, thisstrength

  pulsevel=0.d0

  if (myintime.ge.pulsestart(ipulse)) then
     time=myintime-pulsestart(ipulse)

!!# notice denominator here , half the total pulse time: energy at start of pulse 
!!# is omega2-chirp; at 1/4 pulse omega2-chirp/2; at 3/4 omega2+chirp/2; end of pulse 
!!# omega2+chirp.  Therefore with given value of chirp, will span this range of 
!!# energies over FWHM.

     thisomega2=omega2(ipulse)+chirp(ipulse)*(time-pi/omega(ipulse)/2)/(pi/omega(ipulse)/2)   /2

     thisstrength=pulsestrength(ipulse)*(1d0+ramp(ipulse) * &
          (time-pi/omega(ipulse)/2)/(pi/omega(ipulse)/2))

     if (time.le.pi/omega(ipulse)) then
        pulsevel = thisstrength * sin(time*omega(ipulse))**2 * &
             sin((time-pi/omega(ipulse)/2)*thisomega2 + phaseshift(ipulse))
     endif
  endif

end function pulsevel


!! wtf with longstep...  why did I do it 2* longstep+1?  So for longstep 0, no constant part of
!!  envelope; for longstep 1, constant part is middle 2/3 of pulse.

function longpulselen(myintime, ipulse)
  use pulse_parameters
  use constant_parameters
  use fileptrmod
  implicit none
  integer,intent(in) :: ipulse
  real*8,intent(in) :: myintime
  real*8 :: time, fac, fac2
  real*8 :: longpulselen

  longpulselen=0.d0

  if (chirp(ipulse).ne.0d0.or.ramp(ipulse).ne.0d0) then
     OFLWR "Chirp not supported for length", chirp(ipulse), ipulse; CFLST
  endif
  if (myintime.ge.pulsestart(ipulse)) then
     time=myintime-pulsestart(ipulse)

     if (time.le.pi/omega(ipulse)) then

        if ( (time.le.pi/2.d0/omega(ipulse)/(2*longstep(ipulse)+1)) .or.&
             (time.ge.pi/omega(ipulse) - pi/2.d0/omega(ipulse)/(2*longstep(ipulse)+1)) ) then
           fac=2*omega(ipulse)*(2*longstep(ipulse)+1)*sin(time*omega(ipulse)*(2*longstep(ipulse)+1)) * &
                cos(time*omega(ipulse)*(2*longstep(ipulse)+1))
           fac2=sin(time*omega(ipulse)*(2*longstep(ipulse)+1))**2
        else
           fac=0.d0
           fac2=1.d0
        endif

        longpulselen = pulsestrength(ipulse) * ( &
             fac * sin(time*omega2(ipulse) + phaseshift(ipulse)) &
          + fac2 * omega2(ipulse) * cos(time*omega2(ipulse) + phaseshift(ipulse)) )
     endif
  endif

end function longpulselen


function longpulsevel(myintime, ipulse)
  use pulse_parameters
  use constant_parameters
  implicit none
  integer,intent(in) :: ipulse
  real*8,intent(in) :: myintime
  real*8 :: time, thisomega2
  real*8 :: longpulsevel,thisstrength

  longpulsevel=0.d0

  if (myintime.ge.pulsestart(ipulse)) then
     time=myintime-pulsestart(ipulse)

     if (time.le.pi/omega(ipulse)) then

!! this is right I think  goes to chirp/2  when 2*logstep+1 / 4*longstep+2 -way before half time

        thisomega2=omega2(ipulse)+chirp(ipulse)/2 *(time-pi/omega(ipulse)/2)/&
             (pi/omega(ipulse)/2*(2*longstep(ipulse)+1)/(4*longstep(ipulse)+2))
        thisstrength=pulsestrength(ipulse)*(1d0 +ramp(ipulse) * &
             (time-pi/omega(ipulse)/2)/(pi/omega(ipulse)/2))

!!(2*longstep(ipulse)+1)/(4*longstep(ipulse)+2)))

        if ( (time.le.pi/2.d0/omega(ipulse)/(2*longstep(ipulse)+1)) .or. &
             (time.ge.pi/omega(ipulse) - pi/2.d0/omega(ipulse)/(2*longstep(ipulse)+1)) ) then
          longpulsevel = thisstrength  * &
               sin((time-pi/omega(ipulse)/2)*thisomega2 + phaseshift(ipulse)) * &
               sin(time*omega(ipulse)*(2*longstep(ipulse)+1))**2 
        else
           longpulsevel = thisstrength * &
                sin((time-pi/omega(ipulse)/2)*thisomega2 + phaseshift(ipulse))
        endif
     endif
  endif
  
end function longpulsevel

end module pulsesubmod
