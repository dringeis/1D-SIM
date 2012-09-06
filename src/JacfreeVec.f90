subroutine JacfreeVec (v, Jv, F_uk1, uk1, upts, tauair, epsilon)

  use size
  use global_var
  use option
  
  implicit none
  
  integer :: i
  double precision, intent(in) :: v(1:nx+1), F_uk1(1:nx+1), uk1(1:nx+1)
  double precision, intent(in) :: upts(1:nx+1), tauair(1:nx+1)
  double precision, intent(out):: Jv(1:nx+1)
  double precision, intent(in) :: epsilon
  double precision :: zeta(0:nx+1), eta(0:nx+1)
  double precision :: Cw(1:nx+1), Fpos(1:nx+1)
  double precision :: upos(1:nx+1), b(1:nx+1)

!  double precision xpos(nvar), xneg(nvar), x(nvar),rhs(nvar)

!  double precision Fpos(nvar),Fneg(nvar)
!  double precision epsilon,v(nvar),Jv(nvar)
      
  do i=1, nx+1

     upos(i) = uk1(i) + epsilon * v(i)
!     xneg(i) = x(i) - epsilon * v(i)

  enddo
        
  if (IMEX .eq. 2) then ! IMEX method 2 only
    h=hold
    A=Aold
    call advection (upos, h, A) ! advection scheme for tracers
    call ice_strength () ! Pp_half is Pp/2 where Pp is the ice strength
  endif
  
  call viscouscoefficient (upos, zeta, eta)
  call bvect (tauair, upts, b)
  call Cw_coefficient (upos, Cw)
  call Fu (upos, zeta, eta, Cw, b, Fpos)

  do i=1, nx+1

     Jv(i) = ( Fpos(i)-F_uk1(i) ) / epsilon
            
  enddo

  return
end subroutine JacfreeVec
      



