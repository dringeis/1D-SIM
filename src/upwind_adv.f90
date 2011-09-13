subroutine upwind_adv
  use size
  use resolution
  use global_var
  
  implicit none
  
  integer :: i

  double precision :: deltah(1:nx), flux1h, flux2h
  double precision :: deltaA(1:nx), flux1A, flux2A
  
!------------------------------------------------------------------------
!     compute the deltah in cell i from both sides
!------------------------------------------------------------------------

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

     deltah(i) = Deltat*( flux1h - flux2h ) / Deltax
     deltaA(i) = Deltat*( flux1A - flux2A ) / Deltax

  enddo

!------------------------------------------------------------------------
!     update the tracer values
!     (in a separate do-loop to conserve mass)
!------------------------------------------------------------------------
            
  do i = 1, nx

     h(i) = h(i) + deltah(i)    
     h(i) = max(h(i), 0d0)

     A(i) = A(i) + deltaA(i)    
     A(i) = max(A(i), 0d0)
     A(i) = min(A(i), 1d0)     
     
  enddo
  
  return
end subroutine upwind_adv









