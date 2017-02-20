MODULE shallow_water
use size
! etaw    : sea surface height [m]
! uice : ice velocity               [m/s]

!   h(0)--u(1)--h(1)--u(2)--...u(n)--h(n)--u(n+1)--h(n+1)
!   h(0) and h(n+1) are used for open bcs

  IMPLICIT NONE

  DOUBLE PRECISION, allocatable :: etaw(:), etawn1(:), etawn2(:)
  DOUBLE PRECISION, allocatable :: uwn1(:), uwn2(:)
  
END MODULE shallow_water