
#include "Definitions.INC"


program exact
  use parameters
  implicit none
  integer :: i

  pi=4d0*atan(1d0)
  do i=1,SLN
     nullbuff(i:i)=" "
  enddo

  call getinpfile()
  call getparams()
  call exactalloc()
  call schrodinger()

end program exact

