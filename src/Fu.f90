!****************************************************************************
!     calculates Au-b
!****************************************************************************

subroutine Fu (utp, zeta, eta, Cw, b, Fu_vec)
  use size
  use resolution
  use properties
  use global_var
  use option
  use EVP_const

  implicit none
      
  integer :: i

  double precision, intent(in)  :: utp(1:nx+1)
  double precision, intent(in)  :: zeta(0:nx+1), eta(0:nx+1)
  double precision, intent(in)  :: Cw(1:nx+1)
  double precision, intent(in)  :: b(1:nx+1)

  double precision, intent(out) :: Fu_vec(1:nx+1)
  double precision :: h_at_u

  Fu_vec(1)    = 0d0
  Fu_vec(nx+1) = 0d0

  do i = 2, nx

     Fu_vec(i) = 0.0d0

!------------------------------------------------------------------------
!    rhoice*h*du/dt : tendency term, advection of momentum is neglected
!------------------------------------------------------------------------

     h_at_u = ( h(i) + h(i-1) ) / 2d0

     if (implicit_solv) then
     
        Fu_vec(i) = Fu_vec(i) + ( rho * h_at_u * utp(i) ) / Deltat
     
     elseif (.not. implicit_solv) then

        Fu_vec(i) = Fu_vec(i) + one_or_zero*( rho * h_at_u * utp(i) )/ Deltat

     endif

!------------------------------------------------------------------------
!     Cw*u : water drag term
!------------------------------------------------------------------------
     
     Fu_vec(i) = Fu_vec(i) + Cw(i) * utp(i)
     
!------------------------------------------------------------------------
!     -d ( (zeta+eta) du/dx ) / dx : rheology term
!------------------------------------------------------------------------

     Fu_vec(i) = Fu_vec(i) - &

          (zeta(i)+eta(i)) * (utp(i+1)-utp(i))     / Deltax2 + &
          (zeta(i-1)+eta(i-1)) * (utp(i)-utp(i-1)) / Deltax2
     
!------------------------------------------------------------------------
!     -b : forcing term
!------------------------------------------------------------------------
     
     Fu_vec(i) = Fu_vec(i) - b(i)

!------------------------------------------------------------------------
!     dP/dx term for the EVP...recall it is not included in b
!------------------------------------------------------------------------

     if (.not. implicit_solv) then
        
        Fu_vec(i) = Fu_vec(i) + ( P_half(i) - P_half(i-1) ) / Deltax
        
     endif


  enddo

!  print *, 'L2norm =', sqrt(DOT_PRODUCT(Fu_vec,Fu_vec))

  return
end subroutine Fu






