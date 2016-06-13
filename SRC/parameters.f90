
#include "Definitions.INC"

module constant_parameters
  implicit none
  real*8 :: pi
end module constant_parameters


module fileptrmod
  implicit none
  integer :: mpifileptr=6  !! STANDARD OUT
end module fileptrmod


module tol_parameters
  implicit none
  real*8 :: lntol=1d-4
  real*8 :: invtol=1d-8
end module tol_parameters


module parameters
  use fileptrmod
  use constant_parameters
  implicit none
  integer :: numstates=2
!!$  complex*16 :: stateEnergies(MXST)=0d0
!!$  complex*16 :: couplingmat(MXST,MXST)=0d0
  real*8 :: stateEnergies(MXST)=0d0
  real*8 :: couplingmat(MXST,MXST)=0d0

  real*8 :: expotol=1d-7
  real*8 :: par_timestep=0.05d0
  integer :: numsteps=100
  integer :: notiming=2

  integer :: velflag=0

  integer :: exactflag=0

  character(len=SLN) :: inpfile="Input.Inp"

  character(len=SLN) :: nullbuff

end module parameters



module pulse_parameters           !!      NAMELIST PULSE except for conjgpropflag
  integer :: numpulses=1            !!  number of pulses, enter pulsetype, omega, etc. for each
  integer :: pulsetype(100)=1      !!              !!  Pulsetype=1:  A(t) = pulsestrength * sin(w t)^2,
  real*8  :: omega(100)=1.d0        !!              !!  2:  A(t) = strength * sin(w t)^2 
  real*8 :: omega2(100)=1.d0        !!              !!             * sin(w2 t + phaseshift),
  real*8 :: pulsestart(100)=0d0     !!              !!   
  real*8 :: phaseshift(100)=0d0     !!              !!    pulsestart < t < pulsestart + pi/w; 0 otherwise
  real*8 :: chirp(100)=0d0          !!              !!
  real*8 :: ramp(100)=0d0
  real*8 :: longstep(100)=1d0       !!              !!  Pulsetype 3 available: monochromatic, sinesq start+end
  real*8 :: pulsestrength(100)=.5d0 !!            !!  A_0 = E_0/omega (strength of field)  
  real*8 :: intensity(100)= -1.d0   !!              !! overrides pulse strength.  Intensity, 10^16 W cm^-2 
  real*8 :: pulsetheta(100)=0.d0    !!              !!  angle between polarization and bond axis (radians)
  real*8 :: pulsephi(100)=0.d0      !!              !!  polarization in xy plane
end module pulse_parameters






