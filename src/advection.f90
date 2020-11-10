subroutine advection (un1, utp, hn1in, An1in, hn2in, An2in, hout, Aout)
  use size
  use resolution
  use option
  
  implicit none
  
  integer :: i, k

  logical :: limiter, order2

  double precision, intent(in) :: un1(1:nx+1), utp(1:nx+1)
  double precision, intent(in) :: hn1in(0:nx+1), An1in(0:nx+1)
  double precision, intent(in) :: hn2in(0:nx+1), An2in(0:nx+1)
  double precision, intent(out) :: hout(0:nx+1), Aout(0:nx+1)
  double precision :: hstar(0:nx+1), Astar(0:nx+1)
  double precision :: ustar(1:nx+1)
  double precision :: fluxh(1:nx), fluxA(1:nx)
  double precision :: alpham, um, fmh, fmA ! um=u at mid path, fmh=hdu/dx at mid path
  double precision :: fmhprime, fmAprime
  double precision :: hbef, Abef ! init (before) positions of particles in semilag  
  double precision :: upper, lower

  hout(0) = 0d0    ! closed b.c.s
  hout(nx+1) = 0d0
  Aout(0) = 0d0
  Aout(nx+1) = 0d0
 
  if (adv_scheme .eq. 'upwind') then
  
!------------------------------------------------------------------------
!     compute RHS of dh/dt=-d(hu)/dx (same idea for A)
!------------------------------------------------------------------------

  call fluxh_A (utp, hn1in, An1in, fluxh, fluxA)

!------------------------------------------------------------------------
!     update the tracer values
!     (in a separate do-loop to conserve mass)
!------------------------------------------------------------------------
            
  do i = 1, nx

     hout(i) = hn1in(i) - DtoverDx*fluxh(i)    
     hout(i) = max(hout(i), 0d0)

     Aout(i) = An1in(i) - DtoverDx*fluxA(i)    
     Aout(i) = max(Aout(i), 0d0)
     Aout(i) = min(Aout(i), 1d0)     
     
  enddo
  
  elseif (adv_scheme .eq. 'upwindRK2') then 
  
  call fluxh_A (un1, hn1in, An1in, fluxh, fluxA) 
  
  do i = 1, nx ! predictor step

     hstar(i) = hn1in(i) - (DtoverDx/2d0)*fluxh(i)
     hstar(i) = max(hstar(i), 0d0)

     Astar(i) = An1in(i) - (DtoverDx/2d0)*fluxA(i)    
     Astar(i) = max(Astar(i), 0d0)
     Astar(i) = min(Astar(i), 1d0)     
     
  enddo
  
  ustar = ( utp + un1 ) / 2d0
  call fluxh_A (ustar, hstar, Astar, fluxh, fluxA) 
  
  do i = 1, nx ! corrector step

     hout(i) = hn1in(i) - DtoverDx*fluxh(i)
     hout(i) = max(hout(i), 0d0)

     Aout(i) = An1in(i) - DtoverDx*fluxA(i)    
     Aout(i) = max(Aout(i), 0d0)
     Aout(i) = min(Aout(i), 1d0)     
     
  enddo

  elseif (adv_scheme .eq. 'semilag') then

!------------------------------------------------------------------------ 
!     Semi-Lagrangian scheme for advection. 
!     This is a 3 time level scheme (h is obtained from hn1 and hn2)
!     
!     Staniforth and Côté, Monthly Weather Review 1991.
!     Pellerin et al, Monthly Weather Review 1995. 
!
!------------------------------------------------------------------------ 
     
     limiter=.true. ! see Pellerin et al. MWR 1995
     order2=.true.
     alpham=0.01
     do i = 1, nx

!------------------------------------------------------------------------  
! find velocity at x-alpham  and t=n1
!------------------------------------------------------------------------
     
        do k = 1, 5
           um = (un1(i+1)+un1(i))/2d0 - (un1(i+1)-un1(i))*alpham/Deltax
           if (order2 .and. i .gt. 1 .and. i .lt. nx) then
              um=um + (alpham**2d0)*(un1(i+2)-un1(i+1)-un1(i)+un1(i-1))/(4d0*Deltax2)
           endif
           alpham=Deltat*um
        enddo

!------------------------------------------------------------------------
! find hbef and Abef (initial position of particle at time level n-2)
!------------------------------------------------------------------------
        
        if (i .eq. 1) then
           hbef = hn2in(i) - 2d0 * alpham * ( hn2in(i+1) - hn2in(i) ) / Deltax
           Abef = An2in(i) - 2d0 * alpham * ( An2in(i+1) - An2in(i) ) / Deltax
           if (limiter) then
              upper=max(hn2in(i), hn2in(i+1))
              lower=min(hn2in(i), hn2in(i+1))
              if (hbef .gt. upper) hbef=upper
              if (hbef .lt. lower) hbef=lower
              upper=max(An2in(i), An2in(i+1))
              lower=min(An2in(i), An2in(i+1))
              if (Abef .gt. upper) Abef=upper
              if (Abef .lt. lower) Abef=lower
           endif
        elseif (i .eq. nx) then
           hbef = hn2in(i) - 2d0 * alpham * ( hn2in(i) - hn2in(i-1) ) / Deltax
           Abef = An2in(i) - 2d0 * alpham * ( An2in(i) - An2in(i-1) ) / Deltax
           if (limiter) then
              upper=max(hn2in(i-1), hn2in(i))
              lower=min(hn2in(i-1), hn2in(i))
              if (hbef .gt. upper) hbef=upper
              if (hbef .lt. lower) hbef=lower
              upper=max(An2in(i-1), An2in(i))
              lower=min(An2in(i-1), An2in(i))
              if (Abef .gt. upper) Abef=upper
              if (Abef .lt. lower) Abef=lower
           endif
        else
           hbef = hn2in(i) - ( hn2in(i+1) - hn2in(i-1) )*alpham / Deltax
           Abef = An2in(i) - ( An2in(i+1) - An2in(i-1) )*alpham / Deltax
           
           if (order2) then
              hbef=hbef + 2d0*(alpham**2d0)* &
                              (hn2in(i-1)-2d0*hn2in(i)+hn2in(i+1))/Deltax2
              Abef=Abef + 2d0*(alpham**2d0)* &
                              (An2in(i-1)-2d0*An2in(i)+An2in(i+1))/Deltax2
           endif

           if (limiter) then
              upper=max(hn2in(i-1), hn2in(i), hn2in(i+1))
              lower=min(hn2in(i-1), hn2in(i), hn2in(i+1))
              if (hbef .gt. upper) hbef=upper
              if (hbef .lt. lower) hbef=lower
              upper=max(An2in(i-1), An2in(i), An2in(i+1))
              lower=min(An2in(i-1), An2in(i), An2in(i+1))
              if (Abef .gt. upper) Abef=upper
              if (Abef .lt. lower) Abef=lower
           endif
        endif
        hbef = max(hbef, 0d0)
        Abef = max(Abef, 0d0)
        Abef = min(Abef, 1d0)

!------------------------------------------------------------------------  
! find fmh, fmA (time level n-1)
!------------------------------------------------------------------------ 

        if (i .eq. 1) then
           fmhprime= ( hn1in(i+1)*(un1(i+2)-un1(i+1)) - &
                       hn1in(i)  *(un1(i+1)-un1(i)  ) ) / Deltax2
           fmAprime= ( An1in(i+1)*(un1(i+2)-un1(i+1)) - &
                       An1in(i)  *(un1(i+1)-un1(i)  ) ) / Deltax2
        elseif (i .eq. nx) then
           fmhprime= ( hn1in(i)  *(un1(i+1)-un1(i) ) - &
                       hn1in(i-1)*(un1(i)-un1(i-1) ) ) / Deltax2
           fmAprime= ( An1in(i)  *(un1(i+1)-un1(i) ) - &
                       An1in(i-1)*(un1(i)-un1(i-1) ) ) / Deltax2
        else
           fmhprime=( hn1in(i+1)*(un1(i+2) - un1(i+1)) - &
                      hn1in(i-1)*(un1(i)   - un1(i-1)) ) / (2d0*Deltax2)

           fmAprime=( An1in(i+1)*(un1(i+2) - un1(i+1)) - &
                      An1in(i-1)*(un1(i)   - un1(i-1)) ) / (2d0*Deltax2)
        endif

        fmh = hn1in(i)*(un1(i+1)-un1(i))/Deltax - alpham*fmhprime
        fmA = An1in(i)*(un1(i+1)-un1(i))/Deltax - alpham*fmAprime

        if (order2) then
           
           fmh=fmh+(alpham**2d0)*(hn1in(i+1)*(un1(i+2)-un1(i+1)) - &
                                  2d0*hn1in(i)*(un1(i+1)-un1(i)) + &
                                  hn1in(i-1)*(un1(i)-un1(i-1)) ) / (2d0*(Deltax**3))

           fmA=fmA+(alpham**2d0)*(An1in(i+1)*(un1(i+2)-un1(i+1)) - &
                                  2d0*An1in(i)*(un1(i+1)-un1(i)) + &
                                  An1in(i-1)*(un1(i)-un1(i-1)) ) / (2d0*(Deltax**3))

        endif

!------------------------------------------------------------------------
! find hout, Aout (after, time level n)
!------------------------------------------------------------------------ 

        hout(i) = hbef - 2d0*Deltat*fmh
        hout(i) = max(hout(i), 0d0)
        Aout(i) = Abef - 2d0*Deltat*fmA
        Aout(i) = max(Aout(i), 0d0)
        Aout(i) = min(Aout(i), 1d0)

     enddo


  endif
  
  return
end subroutine advection


subroutine fluxh_A (utp, htp, Atp, fluxh, fluxA)
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

     fluxh(i) = flux2h - flux1h 
     fluxA(i) = flux2A - flux1A 

  enddo
  
end subroutine fluxh_A
  








