subroutine sigma_s(utp, sigma)
  use size
  use resolution
  use numerical
  use rheology
  use EVP_const
  use global_var

  implicit none

  integer :: i
  
  double precision, intent(in) :: utp(1:nx+1)
  double precision, intent(inout):: sigma(0:nx+1)
  double precision :: dudx, deno

  do i = 1, nx ! for tracer points
     
     dudx = ( utp(i+1) - utp(i) ) / Deltax

     deno = alpha*(abs(dudx))

     deno = max( deno, mindeno )
     
     sigma(i) = beta_1 * ( P_half(i)*(1d0/T)*(dudx/deno - 0.5d0) &
                + sigma(i)/Deltate )

      ! check if it should be Pp or P_half...

  enddo

  return
end subroutine sigma_s



