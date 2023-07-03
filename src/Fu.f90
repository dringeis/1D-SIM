!****************************************************************************
!     calculates F = phdu/dt - R
!****************************************************************************

subroutine Fu (utp, un1, un2, htp, R_uk1, Fu_vec)
  use size
  use resolution
  use properties
  use global_var
  use option

  implicit none

  integer :: i

  double precision, intent(in)  :: utp(1:nx+1), un1(1:nx+1), un2(1:nx+1)
  double precision, intent(in)  :: R_uk1(1:nx+1), htp(0:nx+1)

  double precision, intent(out) :: Fu_vec(1:nx+1)
  double precision :: h_at_u

  Fu_vec(1)    = 0d0
  Fu_vec(nx+1) = 0d0

    do i = 2, nx

      Fu_vec(i) = 0.0d0

!------------------------------------------------------------------------
!    rhoice*h*du/dt : tendency term, advection of momentum is neglected
!------------------------------------------------------------------------

      h_at_u = ( htp(i) + htp(i-1) ) / 2d0
      if ( BDF2 .eq. 0 ) then
    Fu_vec(i) = Fu_vec(i) + ( rho * h_at_u * (utp(i)-un1(i)) ) / Deltat
      elseif ( BDF2 .eq. 1 ) then
    Fu_vec(i) = Fu_vec(i) + &
             (rho * h_at_u / (2d0*Deltat)) * ( 3d0*(utp(i)-un1(i)) - (un1(i)-un2(i)) )
      endif

!------------------------------------------------------------------------
!     Substract the R vector
!------------------------------------------------------------------------

      Fu_vec(i) = scaling(i) * ABS( Fu_vec(i) - R_uk1(i) )

    enddo

  return
end subroutine Fu

!****************************************************************************
!     calculates R in du/dt = R/ph
!****************************************************************************

subroutine calc_R (utp, zeta, eta, Cw, Cb, tauair, R_vec)
  use size
  use resolution
  use properties
  use numerical
  use global_var
  use shallow_water
  use option

  implicit none

  integer :: i

  double precision, intent(in)  :: utp(1:nx+1)
  double precision, intent(in)  :: zeta(0:nx+1), eta(0:nx+1)
  double precision, intent(in)  :: Cw(1:nx+1), Cb(1:nx+1), tauair(1:nx+1)

  double precision, intent(out) :: R_vec(1:nx+1)
  double precision :: a_at_u, h_at_u

  R_vec(1)    = 0d0
  R_vec(nx+1) = 0d0

  do i = 2, nx

     R_vec(i) = 0.0d0
     a_at_u = ( A(i) + A(i-1) ) / 2d0
     a_at_u=max(a_at_u, smallA)
     h_at_u = ( h(i) + h(i-1) ) / 2d0

!------------------------------------------------------------------------
!     tauair : air drag term
!------------------------------------------------------------------------

     R_vec(i) = R_vec(i) + a_at_u*tauair(i)

!------------------------------------------------------------------------
!     Cw*u : water drag term
!------------------------------------------------------------------------

!     R_vec(i) = R_vec(i) - a_at_u*Cw(i) * ( utp(i) - uw(i) )
     R_vec(i) = R_vec(i) - Cw(i) * ( utp(i) - uwn2(i) ) ! to be consistent
                                                               ! with NEMO
!------------------------------------------------------------------------
!     Cb*u : bottom drag
!------------------------------------------------------------------------

     ! R_vec(i) = R_vec(i) - Cb(i) * utp(i)

!------------------------------------------------------------------------
!     -rhoh detaw/dx : ocean tilt term
!------------------------------------------------------------------------

     ! R_vec(i) = R_vec(i) - rho * h_at_u * ge * ( etawn1(i) - etawn1(i-1) ) / Deltax

!------------------------------------------------------------------------
!     d ( (zeta+eta) du/dx ) / dx - 1/2dP/dx : rheology term
!------------------------------------------------------------------------

     R_vec(i) = R_vec(i) + &

          (zeta(i)+eta(i)) * (utp(i+1)-utp(i))     / Deltax2 - &
          (zeta(i-1)+eta(i-1)) * (utp(i)-utp(i-1)) / Deltax2

     R_vec(i) = R_vec(i) - ( P_half(i) - P_half(i-1) ) / Deltax

  enddo

  return
end subroutine calc_R



