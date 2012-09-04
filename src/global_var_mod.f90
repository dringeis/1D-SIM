MODULE global_var
use size
! h    : average sea ice thickness  [m]
! uice : ice velocity               [m/s]

!   h(0)--u(1)--h(1)--u(2)--...u(n)--h(n)--u(n+1)--h(n+1)
!   h(0) and h(n+1) are used for open bcs

  IMPLICIT NONE
  DOUBLE PRECISION :: h(0:nx+1)
  DOUBLE PRECISION :: A(0:nx+1)
  DOUBLE PRECISION :: Pp_half(0:nx+1)
  DOUBLE PRECISION :: P_half(0:nx+1)
  DOUBLE PRECISION :: u(1:nx+1)

END MODULE global_var
