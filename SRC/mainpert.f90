
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

