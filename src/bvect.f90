subroutine bvect(tauair, upts, b)
  use size
  use resolution
  use properties
  use global_var
  use option
  use EVP_const

  implicit none
      
  integer :: i
  double precision h_at_u
  double precision, intent(in) :: tauair(1:nx+1)
  double precision, intent(in) :: upts(1:nx+1)
  double precision, intent(out):: b(1:nx+1)

  b(1)    = 0d0 ! close bc
  b(nx+1) = 0d0 ! close bc

  do i = 2, nx

     h_at_u = ( h(i) + h(i-1) ) / 2d0

     if (implicit) then

        b(i) = tauair(i)  - ( p_half(i) - p_half(i-1) ) / Deltax +  &
             ( rho * h_at_u * upts(i) ) / Deltat

     elseif (.not. implicit) then
        
        b(i) = tauair(i) + one_or_zero * ( rho * h_at_u * upts(i) ) / Deltat

     endif
     
  enddo

  return
end subroutine bvect
    
