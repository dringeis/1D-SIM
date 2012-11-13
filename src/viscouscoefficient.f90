subroutine viscouscoefficient(utp, zeta, eta)
  use size
  use resolution
  use rheology
  use option
  use global_var

  implicit none

  integer :: i
  
  double precision, intent(in) :: utp(1:nx+1)
  double precision, intent(out):: zeta(0:nx+1), eta(0:nx+1)
  double precision :: dudx, deno

  do i = 1, nx ! for tracer points
     
     dudx = ( utp(i+1) - utp(i) ) / Deltax

     deno = alpha*sqrt( (dudx)**2d0 + small2 )
!     deno = alpha*(abs(dudx))
!     deno = max( deno, 1d-30 )
     
     zeta(i) = zmax_par * Pp_half(i) * (tanh(1/(zmax_par*deno))) &
          + zetamin

     eta(i)  = zeta(i) * e_2

     if (rep_closure) then  ! replacement closure (Kreysher et al. 2000)
	P_half(i) = zeta(i)*deno 
     endif

  enddo

  if (linear_viscous) then
     zeta = 1d08
     eta  = e_2 * 1d08
  endif

  return
end subroutine viscousCoefficient



