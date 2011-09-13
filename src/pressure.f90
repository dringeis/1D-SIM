subroutine pressure ( )
  use size
  use rheology
  use global_var
  
  implicit none

  integer :: i

  p_half(0)    = 0d0 ! ! sea ice pressure / 2d0
  p_half(nx+1) = 0d0

  do i = 1, nx
     p_half(i) = 0.5d0 * Pstar * h(i) * dexp(-C * ( 1d0 - A(i) ) )
  enddo

!------- set p = 0 at open boundaries for proper care of open bc --------------
!                    see p.1241-1242 for details              
!--- set dh/dx, dA/dx = 0 at the outside cell when there is an open bc --------

  return
end subroutine pressure
      




