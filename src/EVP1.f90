!****************************************************************************
! calculates utp = M^-1(rhs) where rhs is wk1. du=utp= 0 is the initial
! guess when used as a precond. utp is the solution sent back. 
! see p. 3-21 and 3-22
!****************************************************************************

subroutine EVP1 (rhs, utp, zeta, eta, Cw, p_flag, ts)
  use size
  use resolution
  use properties
  use global_var
  use numerical
  use rheology
  use EVP_const

  implicit none
      
  integer :: i, p
  integer, intent(in) :: ts
  logical, intent(in) :: p_flag
  double precision, intent(in)  :: rhs(1:nx+1)
  double precision, intent(in)  :: Cw(1:nx+1)
  double precision, intent(in)  :: zeta(0:nx+1), eta(0:nx+1)
  double precision, intent(out) :: utp(1:nx+1)
  double precision :: sigma(0:nx+1)
  double precision :: h_at_u, B1, gamma, Eyoung, right, left

!------------------------------------------------------------------------
! initial value of sigma and utp
!------------------------------------------------------------------------
  
  if (p_flag) then
     utp = 0d0   ! initial guess for precond (in fact du)
     sigma = 0d0 ! initial guess for precond (in fact dsigma)
  elseif (.not. p_flag) then

     do i = 1, nx ! initial values of sigma

        sigma(i) = (eta(i)+ zeta(i))*( utp(i+1) - utp(i) ) / Deltax &
                  - p_half(i)
     enddo

  endif

  do p = 1, N_sub ! subcycling loop

!------------------------------------------------------------------------
! time step velocity
!------------------------------------------------------------------------

     do i = 2, nx
     
!------------------------------------------------------------------------
!     B1: air drag
!------------------------------------------------------------------------
        
        B1 = rhs(i)

!------------------------------------------------------------------------
!     B1: rho*h*du^p-1 / Deltate
!------------------------------------------------------------------------

        h_at_u = ( h(i) + h(i-1) ) / 2d0

        B1 = B1 + ( rho * h_at_u * utp(i) ) / Deltate
     
!------------------------------------------------------------------------
!     B1: dsigma/dx
!------------------------------------------------------------------------
     
        B1 = B1 + ( sigma(i) - sigma(i-1) )/ Deltax 
     
!------------------------------------------------------------------------
!     advance u from u^p-1 to u^p
!------------------------------------------------------------------------

        gamma = ( rho * h_at_u )*(1d0/Deltate + 1d0/Deltat) + Cw(i)
        
        utp(i) = B1 / gamma

     enddo
     

!------------------------------------------------------------------------
! time step sigma
!------------------------------------------------------------------------
 
     do i = 1, nx ! for tracer points

        Eyoung = Eyoung0 * max(h(i), 0.0001d0)

        right = ( utp(i+1) - utp(i) ) / Deltax + sigma(i) / (Deltate*Eyoung)
        
        left  = 1d0/(Deltate*Eyoung) + ( 1d0 )/( eta(i)+ zeta(i))

        sigma(i) = right / left

     enddo

!     if (ts .gt. 4) then
!        print *, 'truc', p, du(50), dsigma(50)
!     endif

  enddo

  return
end subroutine EVP1






