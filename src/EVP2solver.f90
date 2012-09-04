!****************************************************************************
! calculates utp = M^-1(rhs) where rhs is wk1. du=utp= 0 is the initial
! guess when used as a precond. utp is the solution sent back. 
! see p. 3-21 and 3-22
!****************************************************************************

subroutine EVP2solver (rhs, utp, zeta, eta, Cw, ts)
  use size
  use resolution
  use properties
  use global_var
  use numerical
  use rheology
  use EVP_const

  implicit none
      
  integer :: i, s
  integer, intent(in) :: ts
  double precision, intent(in)  :: rhs(1:nx+1)
  double precision, intent(inout) :: utp(1:nx+1)
  double precision :: Cw(1:nx+1), F_uk1(1:nx+1), us1(1:nx+1)
  double precision :: sigma(0:nx+1), zeta(0:nx+1), eta(0:nx+1), h_at_u(0:nx+1)
  double precision :: B1, gamma, right, left
  double precision :: Fevp(1:nx+1), L2normb ! calc EVP L2norm

  us1 = utp ! at this point us1=us=u^0

  left  = 1d0/(Deltate) + ( 1d0 )/(T*alpha2) ! no change during subcycling

!------------------------------------------------------------------------
! initial value of sigma and utp
!------------------------------------------------------------------------
  
  do i = 1, nx ! initial values of sigma

     sigma(i) = (eta(i)+ zeta(i))*( utp(i+1) - utp(i) ) / Deltax - P_half(i)

  enddo

  do i = 1, nx ! could be improved for precond...

     h_at_u(i) = ( h(i) + h(i-1) ) / 2d0 ! no change during subcycling
      
  enddo

!------------------------------------------------------------------------
! beginning of subcycling loop 
!------------------------------------------------------------------------

  do s = 1, N_sub ! subcycling loop

     if (s .gt. 1) then
        call viscouscoefficient (utp, zeta, eta)
        call Cw_coefficient (utp, Cw)
     endif

!------------------------------------------------------------------------
! calculation of L2norm at each subcycling step
!------------------------------------------------------------------------

!     call Fu (utp, zeta, eta, Cw, rhs, F_uk1)   ! u is u^k-1
     call Fu_EVP (utp, us1, zeta, eta, Cw, rhs, Fevp)   ! u is u^k-1
!     L2norm = sqrt(DOT_PRODUCT(F_uk1,F_uk1))
     L2normb = sqrt(DOT_PRODUCT(Fevp,Fevp))
     print *, 'L2 norm after s subcycles is', ts, s-1, L2normb
!     if (s .eq. 1) nl_target = gamma_nl*L2norm
!     if (L2norm .lt. nl_target) exit

!------------------------------------------------------------------------
! time step sigma
!------------------------------------------------------------------------
 
     do i = 1, nx ! for tracer points

        right = (zeta(i)/ T )*( utp(i+1) - utp(i) ) / Deltax &
              + sigma(i) / Deltate - P_half(i) / (T*alpha2)! could be impr.

        sigma(i) = right / left

     enddo

!------------------------------------------------------------------------
! time step velocity
!------------------------------------------------------------------------

     us1 = utp

     do i = 2, nx
     
!------------------------------------------------------------------------
!     B1: air drag
!------------------------------------------------------------------------
        
        B1 = rhs(i)

!------------------------------------------------------------------------
!     B1: rho*h*du^p-1 / Deltate
!------------------------------------------------------------------------

        B1 = B1 + ( rho * h_at_u(i) * utp(i) ) / Deltate

!------------------------------------------------------------------------
!     B1: rho*h*du^p-1 / Deltat ! to match implicit solution
!------------------------------------------------------------------------

        B1 = B1 - one_or_zero * ( rho * h_at_u(i) * utp(i) ) / Deltat
     
!------------------------------------------------------------------------
!     B1: dsigma/dx
!------------------------------------------------------------------------
     
        B1 = B1 + ( sigma(i) - sigma(i-1) )/ Deltax 
     
!------------------------------------------------------------------------
!     advance u from u^p-1 to u^p
!------------------------------------------------------------------------

        gamma = ( rho * h_at_u(i) )*(1d0/Deltate) + Cw(i)
        
        utp(i) = B1 / gamma

     enddo

  enddo

  return
end subroutine EVP2solver






