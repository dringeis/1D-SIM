subroutine advection (un1, utp, hn1in, An1in, hn2in, An2in, hout, Aout)
  use size
  use resolution
  use option
  
  implicit none
  
  integer :: i, k, lim_scheme, order, ibeg, iend

  logical :: limiter

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
  double precision :: upper, lower, apply_lim1, apply_lim2, cubic_interp, xdist

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
     lim_scheme=2   ! 1: simple, 2: Pellerin et al. MWR 1995
     order=3
     
     if (order .le. 2) then
        ibeg=1
        iend=nx
     elseif (order .eq. 3) then
        ibeg=2
        iend=nx-1
     endif

     alpham=0.01
     do i = 1, nx

!------------------------------------------------------------------------  
! find velocity at x-alpham  and t=n1
!------------------------------------------------------------------------
     
        do k = 1, 5
           um = (un1(i+1)+un1(i))/2d0 - (un1(i+1)-un1(i))*alpham/Deltax
           if (order .gt. 1 .and. i .gt. 1 .and. i .lt. nx) then ! O2...O3 not coded yet  
              um=um + (alpham**2d0)*(un1(i+2)-un1(i+1)-un1(i)+un1(i-1))/(4d0*Deltax2)
           endif
           alpham=Deltat*um
        enddo

!------------------------------------------------------------------------
! find hbef and Abef (initial position of particle at time level n-2)
!------------------------------------------------------------------------
        
        if (i .le. ibeg) then
           hbef = hn2in(i) - 2d0 * alpham * ( hn2in(i+1) - hn2in(i) ) / Deltax
           Abef = An2in(i) - 2d0 * alpham * ( An2in(i+1) - An2in(i) ) / Deltax
           if (limiter) then
              upper=max(hn2in(i), hn2in(i+1))
              lower=min(hn2in(i), hn2in(i+1))
              hbef=apply_lim1(hbef, upper, lower)

              upper=max(An2in(i), An2in(i+1))
              lower=min(An2in(i), An2in(i+1))
              Abef=apply_lim1(Abef, upper, lower)
           endif

         elseif (i .ge. iend) then
           hbef = hn2in(i) - 2d0 * alpham * ( hn2in(i) - hn2in(i-1) ) / Deltax
           Abef = An2in(i) - 2d0 * alpham * ( An2in(i) - An2in(i-1) ) / Deltax
           if (limiter) then
              upper=max(hn2in(i-1), hn2in(i))
              lower=min(hn2in(i-1), hn2in(i))
              hbef=apply_lim1(hbef, upper, lower)

              upper=max(An2in(i-1), An2in(i))
              lower=min(An2in(i-1), An2in(i))
              Abef=apply_lim1(Abef, upper, lower)
           endif
        else

           if (order .eq. 1) then

              hbef = hn2in(i) - ( hn2in(i+1) - hn2in(i-1) )*alpham / Deltax
              Abef = An2in(i) - ( An2in(i+1) - An2in(i-1) )*alpham / Deltax
           
           elseif (order .eq. 2) then
              
              hbef=hn2in(i) - (hn2in(i+1) - hn2in(i-1))*alpham / Deltax + &
                   2d0*(alpham**2d0)*(hn2in(i-1)-2d0*hn2in(i)+hn2in(i+1))/Deltax2

              Abef=An2in(i) - (An2in(i+1) - An2in(i-1))*alpham / Deltax + &
                   2d0*(alpham**2d0)*(An2in(i-1)-2d0*An2in(i)+An2in(i+1))/Deltax2
           
           elseif (order .eq. 3) then

              if (alpham .ge. 0d0) then
                 xdist=(Deltax - 2d0*alpham)/Deltax
                 hbef=cubic_interp (hn2in(i-2), hn2in(i-1), hn2in(i), hn2in(i+1), xdist)
                 Abef=cubic_interp (An2in(i-2), An2in(i-1), An2in(i), An2in(i+1), xdist)
              else ! alpham .lt. 0d0
                 xdist=-2d0*alpham/Deltax
                 hbef=cubic_interp (hn2in(i-1), hn2in(i), hn2in(i+1), hn2in(i+2), xdist)
                 Abef=cubic_interp (An2in(i-1), An2in(i), An2in(i+1), An2in(i+2), xdist)
              endif

           endif

           if (limiter) then ! LIMITER COULD BE IMPROVED FOR O3
             
              ! ---- for h --------
              upper=max(hn2in(i-1), hn2in(i), hn2in(i+1))
              lower=min(hn2in(i-1), hn2in(i), hn2in(i+1))

              if (lim_scheme .eq. 1) then
                 hbef=apply_lim1(hbef, upper, lower)
              elseif (lim_scheme .eq. 2) then
                 hbef=apply_lim2(hbef,upper,lower,alpham,hn2in(i-1),hn2in(i),hn2in(i+1))
              endif
              
              ! ---- for A -------- 
              upper=max(An2in(i-1), An2in(i), An2in(i+1))
              lower=min(An2in(i-1), An2in(i), An2in(i+1))
              
              if (lim_scheme .eq. 1) then
                 Abef=apply_lim1(Abef, upper, lower)
              elseif(lim_scheme .eq. 2) then
                 Abef=apply_lim2(Abef,upper,lower,alpham,An2in(i-1),An2in(i),An2in(i+1))
              endif

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

        if (order .eq. 1) then
           fmh = hn1in(i)*(un1(i+1)-un1(i))/Deltax - alpham*fmhprime
           fmA = An1in(i)*(un1(i+1)-un1(i))/Deltax - alpham*fmAprime

        elseif (order .gt. 1) then ! O2...O3 not coded yet
           
           fmh=hn1in(i)*(un1(i+1)-un1(i))/Deltax - alpham*fmhprime + &
                     (alpham**2d0)*(hn1in(i+1)*(un1(i+2)-un1(i+1)) - &
                                  2d0*hn1in(i)*(un1(i+1)-un1(i)) + &
                                  hn1in(i-1)*(un1(i)-un1(i-1)) ) / (2d0*(Deltax**3))

           fmA=An1in(i)*(un1(i+1)-un1(i))/Deltax - alpham*fmAprime + &
                     (alpham**2d0)*(An1in(i+1)*(un1(i+2)-un1(i+1)) - &
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
  
function apply_lim1(var, upper, lower) result(var_lim)
  double precision, intent(in) :: var, upper, lower ! input
  double precision             :: var_lim ! output
  
  var_lim=var
  if (var .gt. upper) var_lim=upper
  if (var .lt. lower) var_lim=lower

end function

function apply_lim2(var, upper, lower, alpham, var_im1, var_i, var_ip1) result(var_lim)
  use resolution
  
  double precision, intent(in) :: var, upper, lower, alpham, var_im1, var_i, var_ip1 ! input
  double precision             :: var_lim ! output                                    
  double precision             :: slope

  var_lim=var

  if (var .gt. upper .or. var .lt. lower) then
     if (alpham .gt. 0) then ! left of var_i
        slope = ( var_i - var_im1 ) / Deltax
     elseif (alpham .le. 0) then ! right of var_i
        slope = ( var_ip1 - var_i ) / Deltax
     endif
     var_lim = var_i - (slope*alpham)
  endif

end function

function cubic_interp (v1, v2, v3, v4, xdist) result(v_interp)
  use resolution

! see https://www.paulinternet.nl/?page=bicubic 
! dist is between 0 and 1

  double precision, intent(in) :: v1, v2, v3, v4, xdist! input                                 
  double precision             :: v_interp ! output                                                                         
  double precision             :: f1_0 ! 1st derivative of f at x=0
  double precision             :: f1_1 ! 1st derivative of f at x=1 (Dx)
  double precision             :: a, b ! parameters for cubic interpolation

!  f1_0 = ( v3 - v1 ) / (2d0*Deltax)
!  f1_1 = ( v4 - v2 ) / (2d0*Deltax)
  f1_0 = ( v3 - v1 ) / 2d0
  f1_1 = ( v4 - v2 ) / 2d0
  a = 2d0*v2 - 2d0*v3 + f1_0 + f1_1
  b = -3d0*v2 + 3d0*v3 - 2d0*f1_0 -f1_1

  v_interp = a*(xdist**3d0) + b*(xdist**2d0) + f1_0*xdist + v2

end function







