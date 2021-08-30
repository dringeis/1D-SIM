!****************************************************************************
!     solves Au=b with the SOR method or Pdu=rhs (as a precond)
!****************************************************************************

subroutine SOR (b, utp, htp, Atp, zeta, eta, Cw, Cb, p_flag, ts)
  use size
  use resolution
  use properties
  use global_var
  use numerical
  use option

  implicit none

  integer :: i, l
  integer, intent(in) :: ts
  logical, intent(in) :: p_flag ! T: precond, F: standard solver
  double precision, intent(inout) :: utp(1:nx+1) !in: ini guess, out: answer
  double precision, intent(in)  :: zeta(0:nx+1), eta(0:nx+1)
  double precision, intent(in)  :: Cw(1:nx+1), Cb(1:nx+1)
  double precision, intent(in)  :: htp(0:nx+1), Atp(0:nx+1)
  double precision, intent(in)  :: b(1:nx+1)
  double precision              :: D(1:nx+1)

  double precision :: h_at_u, a_at_u, B1, residual ,maxerror

  if (p_flag) then
     utp = 0d0              ! initial guess for precond
     maxiteSOR = iteSOR_pre ! nb of ite for precond
  endif

  do i = 2, nx ! D(i) does not change during the SOR iterations

!------------------------------------------------------------------------
!    rhoice*h*du/dt : tendency term, advection of momentum is neglected
!------------------------------------------------------------------------

     h_at_u = ( htp(i) + htp(i-1) ) / 2d0

     if ( BDF2 .eq. 0 ) then
      D(i) = ( rho * h_at_u ) / Deltat
     elseif ( BDF2 .eq. 1 ) then
      D(i) = ( 3d0 * rho * h_at_u ) / ( 2d0*Deltat )
     endif

!------------------------------------------------------------------------
!     Cw*u : water drag term
!------------------------------------------------------------------------

      a_at_u = ( Atp(i) + Atp(i-1) ) / 2d0
      a_at_u=max(a_at_u, smallA)
      D(i) = D(i) + a_at_u*Cw(i)

!------------------------------------------------------------------------
!     Cb*u : bottom drag term
!------------------------------------------------------------------------

      D(i) = D(i) + Cb(i)

!------------------------------------------------------------------------
!     d ( (zeta+eta) du/dx ) / dx : rheology term
!------------------------------------------------------------------------

      D(i) = D(i) + (zeta(i)+eta(i)+zeta(i-1)+eta(i-1)) / Deltax2
      D(i)=scaling(i)*D(i) ! for JFNK

  enddo

  do l = 1, maxiteSOR
    maxerror = 0d0

    do i = 2, nx

!------------------------------------------------------------------------
!     b : forcing term
!------------------------------------------------------------------------

      B1 = b(i)

!------------------------------------------------------------------------
!     -d ( (zeta+eta) du/dx ) / dx : rheology term
!------------------------------------------------------------------------

      B1 = B1 + scaling(i)*((zeta(i)+eta(i))    *utp(i+1) &
          +  (zeta(i-1)+eta(i-1))*utp(i-1)) / Deltax2

      residual = B1/D(i) - utp(i)
      utp(i) = utp(i) + omega * residual

      if (.not. p_flag) then
        if ( abs( residual ) .gt. maxerror ) then
          maxerror = abs( residual )
        endif
      endif

    enddo

    if (.not. p_flag) then
      if ( maxerror .lt. tol_SOR ) then
        print *, 'nb of SOR ite=', l
        exit
      endif
    endif

  enddo

  return
end subroutine SOR

!****************************************************************************
! forms elements of A and solves Au=b with the SOR method
!****************************************************************************

subroutine SOR_A (b, utp, zeta, eta, Cw, k, ts)
  use size
  use resolution
  use properties
  use global_var
  use numerical

  implicit none

  integer :: i, l
  integer, intent(in) :: k,ts
  double precision, intent(inout) :: utp(1:nx+1)
  double precision, intent(in)  :: zeta(0:nx+1), eta(0:nx+1),Cw(1:nx+1)
  double precision, intent(in)  :: b(1:nx+1)

  double precision :: B1, residual ,maxerror
  double precision :: Aleft(1:nx+1), Adiag(1:nx+1), Aright(1:nx+1)

  call formA (utp, zeta, eta, Cw, ts, k, Aleft, Adiag, Aright)

  do l = 1, maxiteSOR

     maxerror = 0d0

     do i = 2, nx

!------------------------------------------------------------------------
!     b : rhs
!------------------------------------------------------------------------

        B1 = b(i)

!------------------------------------------------------------------------
!     off diag A terms
!------------------------------------------------------------------------

        B1 = B1 - Aleft(i)*utp(i-1) - Aright(i)*utp(i+1)

!------------------------------------------------------------------------
!     get latest u
!------------------------------------------------------------------------

        residual = B1/Adiag(i) - utp(i)
        utp(i) = utp(i) + omega * residual

     if ( abs( residual ) .gt. maxerror ) then
             maxerror = abs( residual )
         endif

     enddo

     if ( maxerror .lt. tol_SOR ) exit

  enddo

  return
end subroutine SOR_A

!****************************************************************************
! forms elements of J and solves Jdu=-F with the SOR method
!****************************************************************************

subroutine SOR_J (utp, Futp, zeta, eta, Cw, upts, tauair, k, ts)
  use size
  use resolution
  use properties
  use global_var
  use numerical

  implicit none

  integer :: i, l
  integer, intent(in) :: k,ts
  double precision, intent(inout) :: utp(1:nx+1)
  double precision, intent(in)  :: zeta(0:nx+1), eta(0:nx+1), upts(1:nx+1)
  double precision, intent(in)  :: tauair(1:nx+1), Cw(1:nx+1)
  double precision, intent(in)  :: Futp(1:nx+1)
  double precision              :: du(1:nx+1)

  double precision :: B1, residual ,maxerror
  double precision :: Jleft(1:nx+1), J(1:nx+1), Jright(1:nx+1)

  call formJacobian (utp, Futp, upts, tauair, ts, k, Jleft, J, Jright)

  du = 0d0

  do l = 1, maxiteSOR

     maxerror = 0d0

     do i = 2, nx

!------------------------------------------------------------------------
!     -F : rhs
!------------------------------------------------------------------------

        B1 = -Futp(i)

!------------------------------------------------------------------------
!     off diag Jac terms
!------------------------------------------------------------------------

        B1 = B1 - Jleft(i)*du(i-1) - Jright(i)*du(i+1)

!------------------------------------------------------------------------
!     get latest du
!------------------------------------------------------------------------

        residual = B1/J(i) - du(i)
        du(i) = du(i) + omega * residual

     if ( abs( residual ) .gt. maxerror ) then
             maxerror = abs( residual )
         endif

     enddo

     print *, 'max error', l, maxerror

     if ( maxerror .lt. tol_SOR ) exit

  enddo

! call output_u_and_du ( ts, k, uk1, du )
  utp = utp + du

  return
end subroutine SOR_J

