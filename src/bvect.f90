subroutine bvect(tauair, un1, b)
  use size
  use resolution
  use properties
  use global_var

  implicit none
      
  integer :: i
  double precision h_at_u
  double precision, intent(in) :: tauair(1:nx+1), un1(1:nx+1)
  double precision, intent(out):: b(1:nx+1)

  b(1)    = 0d0 ! close bc
  b(nx+1) = 0d0 ! close bc

  do i = 2, nx

     h_at_u = ( h(i) + h(i-1) ) / 2d0

     b(i) = tauair(i)  - ( P_half(i) - P_half(i-1) ) / Deltax +  &
             ( rho * h_at_u * un1(i) ) / Deltat

  enddo
   
  return
end subroutine bvect
    
