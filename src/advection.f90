subroutine advection (utp, htp, Atp)
  use size
  use resolution
  use option
  
  implicit none
  
  integer :: i

  double precision, intent(in) :: utp(1:nx+1)
  double precision, intent(inout) :: htp(0:nx+1), Atp(0:nx+1)
  double precision :: fluxh(1:nx), flux1h, flux2h
  double precision :: fluxA(1:nx), flux1A, flux2A
 
  htp(0) = 0d0    ! closed b.c.s
  htp(nx+1) = 0d0
  Atp(0) = 0d0
  Atp(nx+1) = 0d0
 
  if (adv_scheme .eq. 'upwind') then
  
!------------------------------------------------------------------------
!     compute RHS of dh/dt=-d(hu)/dx (same idea for A)
!------------------------------------------------------------------------

  call fluxh_A (fluxh, fluxA, utp, htp, Atp)

!------------------------------------------------------------------------
!     update the tracer values
!     (in a separate do-loop to conserve mass)
!------------------------------------------------------------------------
            
  do i = 1, nx

     htp(i) = htp(i) + Deltat*fluxh(i)    
     htp(i) = max(htp(i), 0d0)

     Atp(i) = Atp(i) + Deltat*fluxA(i)    
     Atp(i) = max(Atp(i), 0d0)
     Atp(i) = min(Atp(i), 1d0)     
     
  enddo
  
  elseif (adv_scheme .eq. 'upwindRK2') then 
  stop
  endif
  
  return
end subroutine advection


subroutine fluxh_A (fluxh, fluxA, utp, htp, Atp)
  use size
  use resolution
  
  implicit none
  
  integer :: i

  double precision, intent(in) :: utp(1:nx+1)
  double precision, intent(in) :: htp(0:nx+1), Atp(0:nx+1)
  double precision, intent (out) :: fluxh(1:nx), fluxA(1:nx)
  double precision :: flux1h, flux2h, flux1A, flux2A
  
  do i = 1, nx

! verify CFL condition

     if (utp(i) .gt. Deltax/Deltat) print *, 'WARNING: u > dx/dt', i,utp(i)

! calculate deltah at tracer point i

     if (utp(i) .ge. 0d0) then ! left side of cell
        flux1h = utp(i)*htp(i-1) 
        flux1A = utp(i)*Atp(i-1) 
     else
        flux1h = utp(i)*htp(i)
        flux1A = utp(i)*Atp(i)
     endif
     
     if (utp(i+1) .ge. 0d0) then ! right side of cell
        flux2h = utp(i+1)*htp(i) 
        flux2A = utp(i+1)*Atp(i)
     else
        flux2h = utp(i+1)*htp(i+1)
        flux2A = utp(i+1)*Atp(i+1)
     endif

     fluxh(i) = ( flux1h - flux2h ) / Deltax
     fluxA(i) = ( flux1A - flux2A ) / Deltax

  enddo
  
end subroutine fluxh_A
  







