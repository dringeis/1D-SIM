!****************************************************************************
! EVP solver in 1D. Author JF Lemieux 8 jan 2014.
!****************************************************************************

subroutine EVP2solver (tauair, utp, ts, solver)
  use size
  use resolution
  use properties
  use global_var
  use shallow_water
  use numerical
  use rheology

  implicit none

  integer :: i, s
  integer, intent(in) :: ts, solver
  double precision, intent(in)  :: tauair(1:nx+1)
  double precision, intent(inout) :: utp(1:nx+1)
  double precision :: Cw(1:nx+1), Cb(1:nx+1), F_uk1(1:nx+1), R_uk1(1:nx+1), un1tp(1:nx+1)
  double precision :: zeta(0:nx+1), eta(0:nx+1), h_at_u(0:nx+1), a_at_u(0:nx+1)
  double precision :: B1, gamma, right, left, L2norm
  double precision, save :: sigma(0:nx+1)
!  double precision :: Fevp(1:nx+1), L2normb ! calc EVP L2norm

  left  = 1d0/(Deltate) + ( 1d0 )/(T*alpha2) ! no change during subcycling

!------------------------------------------------------------------------
! initial value of sigma and utp
!------------------------------------------------------------------------

  un1tp = utp ! for EVP* solver

  if ( ts .eq. 1 ) then ! initial sigma set to VP for 1st time level (only) for standard evp

    call viscouscoefficient (utp, zeta, eta)
    call Cw_coefficient (utp, Cw, Cb)
    do i = 1, nx

      sigma(i) = (eta(i)+ zeta(i))*( utp(i+1) - utp(i) ) / Deltax - P_half(i)

    enddo
  endif

  do i = 1, nx ! could be improved for precond...

     h_at_u(i) = ( h(i) + h(i-1) ) / 2d0 ! no change during subcycling
     a_at_u(i) = ( A(i) + A(i-1) ) / 2d0
     a_at_u(i)=max(a_at_u(i), smallA)

  enddo

!------------------------------------------------------------------------
! beginning of subcycling loop
!------------------------------------------------------------------------

  do s = 1, N_sub ! subcycling loop

    call viscouscoefficient (utp, zeta, eta)
    call Cw_coefficient (utp, Cw, Cb)

!------------------------------------------------------------------------
! calculation of L2norm at each subcycling step for EVP* solver
!------------------------------------------------------------------------

    if ( solver .eq. 4 ) then
      call calc_R (utp, zeta, eta, Cw, Cb, tauair, R_uk1)
      call Fu (utp, un1tp, un1tp, h, R_uk1, F_uk1)
      L2norm = sqrt(DOT_PRODUCT(F_uk1,F_uk1))
      print *, 'L2 norm after s subcycles is', ts, s-1, L2norm
     endif

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

     do i = 2, nx

!------------------------------------------------------------------------
!     B1: air drag
!------------------------------------------------------------------------

        B1 = a_at_u(i) * tauair(i)

!------------------------------------------------------------------------
!     B1: part of water drag
!------------------------------------------------------------------------

        B1 = B1 + a_at_u(i) * Cw(i) * uwn2(i) ! to be consistent with NEMO

!------------------------------------------------------------------------
!     -rhoh detaw/dx : ocean tilt term
!------------------------------------------------------------------------

        B1 = B1 - rho * h_at_u(i) * ge * ( etawn1(i) - etawn1(i-1) ) / Deltax

!------------------------------------------------------------------------
!     B1: rho*h*du^p-1 / Deltate
!------------------------------------------------------------------------

        B1 = B1 + ( rho * h_at_u(i) * utp(i) ) / Deltate

!------------------------------------------------------------------------
!     B1: rho*h*u^n-1 / Deltat ! for EVP* solver
!------------------------------------------------------------------------
    if ( solver .eq. 4 ) then
      B1 = B1 + ( rho * h_at_u(i) * un1tp(i) ) / Deltat
    endif
!------------------------------------------------------------------------
!     B1: dsigma/dx
!------------------------------------------------------------------------

        B1 = B1 + ( sigma(i) - sigma(i-1) )/ Deltax

!------------------------------------------------------------------------
!     advance u from u^s-1 to u^s
!------------------------------------------------------------------------
    if ( solver .eq. 3 ) then
      gamma = ( rho * h_at_u(i) )*(1d0/Deltate) + a_at_u(i) * Cw(i) + Cb(i)
    elseif (solver .eq. 4 ) then
      gamma = ( rho * h_at_u(i) )*(1d0/Deltate) + &
               ( rho * h_at_u(i) )*(1d0/Deltat) + a_at_u(i) * Cw(i) + Cb(i)
    endif

      utp(i) = B1 / gamma

    enddo

  enddo

  return
end subroutine EVP2solver






