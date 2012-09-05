subroutine advection
  use size
  use resolution
  use global_var
  use option
  
  implicit none
  
  integer :: i

  double precision :: fluxh(1:nx), flux1h, flux2h
  double precision :: fluxA(1:nx), flux1A, flux2A
  
  if (adv_scheme .eq. 'upwind') then
  
!------------------------------------------------------------------------
!     compute RHS of dh/dt=-d(hu)/dx (same idea for A)
!------------------------------------------------------------------------

  call fluxh_A (fluxh, fluxA)

!------------------------------------------------------------------------
!     update the tracer values
!     (in a separate do-loop to conserve mass)
!------------------------------------------------------------------------
            
  do i = 1, nx

     h(i) = h(i) + Deltat*fluxh(i)    
     h(i) = max(h(i), 0d0)

     A(i) = A(i) + Deltat*fluxA(i)    
     A(i) = max(A(i), 0d0)
     A(i) = min(A(i), 1d0)     
     
  enddo
  
  elseif (adv_scheme .eq. 'upwindRK2') then 
  stop
  endif
  
  return
end subroutine advection


subroutine fluxh_A (fluxh, fluxA)
  use size
  use resolution
  use global_var
  
  implicit none
  
  integer :: i

  double precision, intent (out) :: fluxh(1:nx), fluxA(1:nx)
  double precision :: flux1h, flux2h, flux1A, flux2A
  
  h(0) = 0d0
  h(nx+1) = 0d0
  A(0) = 0d0
  A(nx+1) = 0d0

  do i = 1, nx

! verify CFL condition

     if (u(i) .gt. Deltax/Deltat) print *, 'WARNING: u > dx/dt', i,u(i)

! calculate deltah at tracer point i

     if (u(i) .ge. 0d0) then ! left side of cell
        flux1h = u(i)*h(i-1) 
        flux1A = u(i)*A(i-1) 
     else
        flux1h = u(i)*h(i)
        flux1A = u(i)*A(i)
     endif
     
     if (u(i+1) .ge. 0d0) then ! right side of cell
        flux2h = u(i+1)*h(i) 
        flux2A = u(i+1)*A(i)
     else
        flux2h = u(i+1)*h(i+1)
        flux2A = u(i+1)*A(i+1)
     endif

     fluxh(i) = ( flux1h - flux2h ) / Deltax
     fluxA(i) = ( flux1A - flux2A ) / Deltax

  enddo
  
end subroutine fluxh_A
  








