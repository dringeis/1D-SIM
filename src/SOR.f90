!****************************************************************************
!     solves Au=b with the SOR method or Pdu=rhs (as a precond)
!****************************************************************************

subroutine SOR (rhs, utp, zeta, eta, Cw, p_flag, ts)
  use size
  use resolution
  use properties
  use global_var
  use numerical

  implicit none
      
  integer :: i, l
  integer, intent(in) :: ts
  logical, intent(in) :: p_flag ! T: precond, F: standard solver
  double precision, intent(inout) :: utp(1:nx+1) !in: ini guess, out: answer
  double precision, intent(in)  :: zeta(0:nx+1), eta(0:nx+1)
  double precision, intent(in)  :: Cw(1:nx+1)
  double precision, intent(in)  :: rhs(1:nx+1)
  double precision              :: D(1:nx+1)

  double precision :: h_at_u, B1, residual ,maxerror

  if (p_flag) then
     utp = 0d0              ! initial guess for precond
     maxiteSOR = iteSOR_pre ! nb of ite for precond
  endif

  do i = 2, nx ! D(i) does not change during the SOR iterations

!------------------------------------------------------------------------
!    rhoice*h*du/dt : tendency term, advection of momentum is neglected
!------------------------------------------------------------------------

     h_at_u = ( h(i) + h(i-1) ) / 2d0
     D(i) = ( rho * h_at_u ) / Deltat 

!------------------------------------------------------------------------
!     Cw*u : water drag term
!------------------------------------------------------------------------
     
     D(i) = D(i) + Cw(i)

!------------------------------------------------------------------------
!     d ( (zeta+eta) du/dx ) / dx : rheology term
!------------------------------------------------------------------------

     D(i) = D(i) + (zeta(i)+eta(i)+zeta(i-1)+eta(i-1)) / Deltax2

  enddo

  do l = 1, maxiteSOR
     
     maxerror = 0d0

     do i = 2, nx
        
!------------------------------------------------------------------------
!     b : forcing term
!------------------------------------------------------------------------
        
        B1 = rhs(i)

!------------------------------------------------------------------------
!     -d ( (zeta+eta) du/dx ) / dx : rheology term
!------------------------------------------------------------------------

        B1 = B1 + ((zeta(i)+eta(i))    *utp(i+1) &
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
        if ( maxerror .lt. tol_SOR ) exit
     endif

  enddo

  return
end subroutine SOR






