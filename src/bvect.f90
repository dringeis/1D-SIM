subroutine bvect(tauair, un1, Cw, b)
  use size
  use resolution
  use properties
  use numerical
  use global_var

  implicit none
      
  integer :: i
  double precision h_at_u, a_at_u, ge
  double precision, intent(in) :: tauair(1:nx+1), un1(1:nx+1), Cw(1:nx+1)
  double precision, intent(out):: b(1:nx+1)

  b(1)    = 0d0 ! close bc
  b(nx+1) = 0d0 ! close bc
  ge = 9.8d0
  
  
  do i = 2, nx

     h_at_u = ( h(i) + h(i-1) ) / 2d0
     a_at_u = ( A(i) + A(i-1) ) / 2d0
     a_at_u=max(a_at_u, smallA)

     b(i) = a_at_u*tauair(i) + a_at_u*Cw(i)*uw(i) - &
            ( P_half(i) - P_half(i-1) ) / Deltax + ( rho * h_at_u * un1(i) ) / Deltat - &
            rho * h_at_u * ge * ( etaw(i) - etaw(i-1) ) / Deltax

  enddo
   
  return
end subroutine bvect
    
