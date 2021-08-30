
subroutine calc_scaling ( Atp )

  use size
  use numerical
  use global_var

!--------------------------------------------------------------------------
! Preconditioner: calculates M x = rhs where rhs is wk1. x is the initial
! guess. x is set to 0 here. The solution x is then put in wk2.
!--------------------------------------------------------------------------

  implicit none

  integer :: i
  double precision, intent(in)  :: Atp(0:nx+1)
  double precision :: a_at_u

  do i = 2, nx

     a_at_u = ( Atp(i) + Atp(i-1) ) / 2d0
     a_at_u=max(a_at_u, smallA)
     scaling(i)=1d0/a_at_u

  enddo


end subroutine calc_scaling






