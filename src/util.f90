subroutine meantracer(var, meanvalue)
  use size
  
  implicit none
  
  integer :: i
  double precision, intent(in) :: var(0:nx+1)
  double precision, intent(out) :: meanvalue
  
  meanvalue = 0d0
  
  do i = 1, nx

     meanvalue = meanvalue + var(i)

  enddo

  meanvalue = meanvalue / (nx*1d0)
  
  return
end subroutine meantracer

subroutine check_neg_vel(utp)
  use size

  implicit none

  integer :: i
  double precision :: eps, umin,umax
  double precision, intent(in) :: utp(1:nx+1)

  eps=1d-08
  umin=1d08
  umax=-1d08

  do i = 2, nx ! recall u(1)=0

    if (utp(i) .gt. umax) umax=utp(i)
    if (utp(i) .lt. umin) umin=utp(i)

    if (utp(i) .lt. -eps) then
!       print *, utp(i),'NEGATIVE VELOCITIES HERE AT i=', i
!       stop
    endif

  enddo

  print *, 'min, max u =', umin, umax

  return
end subroutine check_neg_vel

subroutine minmaxtracer(var,id, ts)
  use size

  implicit none

  integer :: i
  integer, save :: negice, negwater
  integer, intent(in) :: id, ts ! h: id=1, A: id=2, zeta: id=3
  double precision, intent(in) :: var(0:nx+1)
  double precision :: vmin, vmax
  
  if (ts .eq. 1) then
    negice=0
    negwater=0
  endif
  
  vmin=1d100
  vmax=-1d100

  do i = 1, nx

     if (var(i) .gt. vmax) vmax=var(i)
     if (var(i) .lt. vmin) vmin=var(i)

  enddo

  if (id .eq. 1) then
     print *, 'min, max h =', vmin, vmax
  elseif (id .eq. 2) then
     print *, 'min, max A =', vmin, vmax
  elseif (id .eq. 3) then
  if (vmin .lt. 0d0) negice=negice+1
     print *, 'min, max u ice =', vmin, vmax, negice
  elseif (id .eq. 4) then
     print *, 'min, max eta water =', vmin, vmax
  elseif (id .eq. 5) then
  if (vmin .lt. 0d0) negwater=negwater+1
     print *, 'min, max u water =', vmin, vmax, negwater
  else
     print *, 'WRONG ID'
  endif

  return
end subroutine minmaxtracer

subroutine stab_condition(Cw,zeta)
  use size
  use global_var
  use rheology
  use resolution
  use properties

  implicit none

  integer :: i, imin

  double precision, intent(in) :: Cw(1:nx+1), zeta(0:nx+1)
  double precision :: LHS1, LHS2, LHS3, LHS4, LHS, LHSmin
  double precision :: m2,mu0, c0, gamma0, nu0 
  double precision :: A_at_u, h_at_u, zeta_at_u
  double precision :: LHStp(4), Amin, hmin

  m2=((2d0*3.1415926)/Deltax)**2d0
  LHSmin=1d100

  do i = 2, nx

     h_at_u = ( h(i-1) + h(i) )/2d0
     h_at_u = max(h_at_u, 1d-09)
     A_at_u = ( A(i-1) + A(i) )/2d0
     zeta_at_u = ( zeta(i-1) + zeta(i) )/2d0

     mu0=0.5d0*Pstar*exp(-C*(1-A_at_u))
     c0=sqrt(mu0*(1d0+C*A_at_u)/rho)
     nu0=alpha2*zeta_at_u/(rho*h_at_u)
     gamma0=Cw(i)/(rho*h_at_u) + m2*nu0

     LHS1=-0.25d0*(m2**2d0)*(c0**4d0)*(Deltat**4d0)
     LHS2=0.5d0*m2*(c0**2d0)*gamma0*(Deltat**3d0)
     LHS3=(m2*(c0**2d0)+0.75d0*gamma0**2d0)*Deltat**2d0
     LHS4=gamma0*Deltat

     LHS=LHS1+LHS2+LHS3+LHS4
     
     if (LHS .lt. LHSmin) then
        LHSmin = LHS
        LHStp(1)=LHS1
        LHStp(2)=LHS2
        LHStp(3)=LHS3
        LHStp(4)=LHS4
        imin=i
        Amin=A_at_u
        hmin=h_at_u
     endif

     if (LHS .lt. 0) then
        print *, 'Dude we are in trouble!!!', Deltat, LHS, LHSmin
      !  stop
     endif

  enddo

  print *, 'ici le min', LHSmin, LHStp(1),LHStp(2),LHStp(3),LHStp(4)

  return
end subroutine stab_condition

!****************************************************************************
!     diagnostic for ice-ocean stress
!****************************************************************************

MODULE diag_stress

IMPLICIT NONE

  DOUBLE PRECISION :: tauaidiag, tauiwdiag, tauwidiag
  
END MODULE diag_stress

subroutine calc_diag_stress (idiag, utp, Cw, tauair)
  use size
  use numerical
  use global_var
  use shallow_water
  use diag_stress

  implicit none
      
  integer, intent(in) :: idiag
  
  double precision, intent(in)  :: utp(1:nx+1), Cw(1:nx+1), tauair(1:nx+1)
  double precision :: a_at_u, h_at_u

  a_at_u = ( A(idiag) + A(idiag-1) ) / 2d0
  a_at_u=max(a_at_u, smallA)

!------------------------------------------------------------------------
!     air stress on ice
!------------------------------------------------------------------------  
  
  tauaidiag = a_at_u*tauair(idiag)
  
!------------------------------------------------------------------------
!     water stress on ice
!------------------------------------------------------------------------
     
  tauwidiag = -1d0 * a_at_u*Cw(idiag) * ( utp(idiag) - uwn2(idiag) ) ! to be consistent
                                                                    ! with NEMO
  return
end subroutine calc_diag_stress








