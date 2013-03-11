!****************************************************************************
!     calculates F = phdu/dt - R 
!****************************************************************************

subroutine Fu (utp, upts, Rpts, R_uk1, Fu_vec)
  use size
  use resolution
  use properties
  use global_var
  use option

  implicit none
      
  integer :: i

  double precision, intent(in)  :: utp(1:nx+1), upts(1:nx+1)
  double precision, intent(in)  :: Rpts(1:nx+1), R_uk1(1:nx+1)

  double precision, intent(out) :: Fu_vec(1:nx+1)
  double precision :: h_at_u, h_at_u_pts

  Fu_vec(1)    = 0d0
  Fu_vec(nx+1) = 0d0

  if ( CN .eq. 0 ) then

    do i = 2, nx

      Fu_vec(i) = 0.0d0

!------------------------------------------------------------------------
!    rhoice*h*du/dt : tendency term, advection of momentum is neglected
!------------------------------------------------------------------------

      h_at_u = ( h(i) + h(i-1) ) / 2d0
      Fu_vec(i) = Fu_vec(i) + ( rho * h_at_u * (utp(i)-upts(i)) ) / Deltat

!------------------------------------------------------------------------
!     Substract the R vector
!------------------------------------------------------------------------
     
      Fu_vec(i) = Fu_vec(i) - R_uk1(i)

    enddo

  elseif ( CN .eq. 1 ) then

    do i = 2, nx

      Fu_vec(i) = 0.0d0

!------------------------------------------------------------------------
!    2*rhoice*h*du/dt : tendency term, advection of momentum is neglected
!
!    2 is hard coded here....watchout (voir p. 42 EC-4)
!------------------------------------------------------------------------

     h_at_u = ( h(i) + h(i-1) ) / 2d0
     h_at_u_pts = ( hpts(i) + hpts(i-1) ) / 2d0
     h_at_u_pts = max ( h_at_u_pts, 1d-25 ) !to avoid div by 0

     Fu_vec(i) = Fu_vec(i) + 2d0*( rho * h_at_u * (utp(i)-upts(i)) ) / Deltat
     
!------------------------------------------------------------------------
!     Substract the R vectors
!------------------------------------------------------------------------
     
     Fu_vec(i) = Fu_vec(i) - R_uk1(i) - h_at_u * Rpts(i) / h_at_u_pts

  enddo
  endif

  return
end subroutine Fu

!****************************************************************************
!     calculates R in du/dt = R/ph
!****************************************************************************

subroutine calc_R (utp, zeta, eta, Cw, tauair, R_vec)
  use size
  use resolution
  use properties
  use global_var
  use option

  implicit none
      
  integer :: i

  double precision, intent(in)  :: utp(1:nx+1)
  double precision, intent(in)  :: zeta(0:nx+1), eta(0:nx+1)
  double precision, intent(in)  :: Cw(1:nx+1), tauair(1:nx+1)

  double precision, intent(out) :: R_vec(1:nx+1)
  
  R_vec(1)    = 0d0
  R_vec(nx+1) = 0d0

  do i = 2, nx

     R_vec(i) = 0.0d0

!------------------------------------------------------------------------
!     tauair : air drag term
!------------------------------------------------------------------------

     R_vec(i) = R_vec(i) + tauair(i)

!------------------------------------------------------------------------
!     Cw*u : water drag term
!------------------------------------------------------------------------
     
     R_vec(i) = R_vec(i) - Cw(i) * utp(i)
     
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



