MODULE shallow_water
use size
! etaw    : sea surface height [m]
! uice : ice velocity               [m/s]

!   h(0)--u(1)--h(1)--u(2)--...u(n)--h(n)--u(n+1)--h(n+1)
!   h(0) and h(n+1) are used for open bcs

  IMPLICIT NONE

  DOUBLE PRECISION, allocatable :: etaw(:), etawn1(:), etawn2(:)
  DOUBLE PRECISION, allocatable :: uwn1(:), uwn2(:)
  DOUBLE PRECISION :: Hw
  
!------------------------------------------------------------------------
!     Subroutine for advecting etaw
!------------------------------------------------------------------------
  
END MODULE shallow_water

subroutine advect_etaw (uwtp)
  
  use size
  use resolution
  use shallow_water
  
  implicit none
  
  double precision, intent(in) :: uwtp(1:nx+1)
  double precision :: RHS, Htleft, Htright ! Ht=Hw+etaw
  integer :: i
  
  ! leapfrog approach with centered (n1) term on the RHS
  
  do i = 1, nx
  
   Htleft  = Hw + ( etawn1(i-1) + etawn1(i) ) /2d0
   Htright = Hw + ( etawn1(i+1) + etawn1(i) ) /2d0
   RHS     = -1d0*( uwtp(i+1)*Htright - uwtp(i)*Htleft ) / Deltax
   etaw(i) = etawn2(i) + 2d0*Deltat*RHS
  
  enddo
  
end subroutine advect_etaw