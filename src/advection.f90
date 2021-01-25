subroutine advection (un1, utp, hn1in, An1in, hn2in, An2in, hout, Aout)
  use size
  use resolution
  use option
  
  implicit none
  
  integer :: i, k, lim_scheme, order, ie, iw, caseSL, CIoption

  logical :: SLlimiter

  double precision, intent(in) :: un1(1:nx+1), utp(1:nx+1)
  double precision, intent(in) :: hn1in(0:nx+1), An1in(0:nx+1)
  double precision, intent(in) :: hn2in(0:nx+1), An2in(0:nx+1)
  double precision, intent(out) :: hout(0:nx+1), Aout(0:nx+1)
  double precision :: hstar(0:nx+1), Astar(0:nx+1)
  double precision :: ustar(1:nx+1)
  double precision :: fluxh(1:nx), fluxA(1:nx), flux, div(nx)
  double precision :: alpham, rhsh, rhsA ! rhsh=hdu/dx at mid path
  double precision :: rhshprime, rhsAprime
  double precision :: fw, fe, fxw, fxe, deno
  double precision :: hbef, Abef ! init (before) positions of particles in semilag  
  double precision :: upper, lower, xd, xdn1, xdn2, uinterp
  double precision :: apply_lim1, apply_lim2, fx, cubic_interp, calc_flux ! functions

  hout(0) = 0d0    ! closed b.c.s
  hout(nx+1) = 0d0
  Aout(0) = 0d0
  Aout(nx+1) = 0d0
 
  if (adv_scheme .eq. 'upwind') then
  
!------------------------------------------------------------------------
!     update h following dh/dt=-d(hu)/dx (same idea for A)
!     (in a separate do-loop to conserve mass)
!------------------------------------------------------------------------
            
     do i = 1, nx

        flux=calc_flux(utp(i),utp(i+1),hn1in(i-1),hn1in(i), hn1in(i+1)) ! for h                                                                
        hout(i) = hn1in(i) - DtoverDx*flux
        hout(i) = max(hout(i), 0d0)

        flux=calc_flux(utp(i),utp(i+1),An1in(i-1),An1in(i), An1in(i+1)) ! for A
        Aout(i) = An1in(i) - DtoverDx*flux
        Aout(i) = max(Aout(i), 0d0)
        Aout(i) = min(Aout(i), 1d0)     
     
     enddo
  
  elseif (adv_scheme .eq. 'upwindRK2') then 
  
     do i = 1, nx ! predictor step

        flux=calc_flux(un1(i),un1(i+1),hn1in(i-1),hn1in(i),hn1in(i+1)) ! for h
        hstar(i) = hn1in(i) - (DtoverDx/2d0)*flux
        hstar(i) = max(hstar(i), 0d0)

        flux=calc_flux(un1(i),un1(i+1),An1in(i-1),An1in(i),An1in(i+1)) ! for A
        Astar(i) = An1in(i) - (DtoverDx/2d0)*flux
        Astar(i) = max(Astar(i), 0d0)
        Astar(i) = min(Astar(i), 1d0)     
     
     enddo
  
     ustar = ( utp + un1 ) / 2d0
  
     do i = 1, nx ! corrector step

        flux=calc_flux(ustar(i),ustar(i+1),hstar(i-1),hstar(i),hstar(i+1)) ! for h  
        hout(i) = hn1in(i) - DtoverDx*flux
        hout(i) = max(hout(i), 0d0)
        
        flux=calc_flux(ustar(i),ustar(i+1),Astar(i-1),Astar(i),Astar(i+1)) ! for A
        Aout(i) = An1in(i) - DtoverDx*flux
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
     
     SLlimiter=.true. ! see Pellerin et al. MWR 1995
     CIoption=2 ! defines region for cubic interpolation of u at mid point
     lim_scheme=2   ! 1: simple, 2: Pellerin et al. MWR 1995
     order=4        ! cubic interp
     
     call calc_div(un1,div) ! calc divergence at n-1 for RHS [hdiv(u)]^{n-1}

     alpham=0.01
     do i = 1, nx
        
        caseSL=1
        if (i .lt. 3 .or. i .gt. nx-2) caseSL=2 ! upwind close to walls

        if (caseSL == 1) then

!------------------------------------------------------------------------  
! find velocity at x-alpham  and t=n1
!------------------------------------------------------------------------

           do k = 1, 5
              if (order == 1) then
                 uinterp = (un1(i+1)+un1(i))/2d0 - (un1(i+1)-un1(i))*alpham/Deltax
              elseif (order == 2) then
                 uinterp = (un1(i+1)+un1(i))/2d0 - (un1(i+1)-un1(i))*alpham/Deltax + &
                      (alpham**2d0)*(un1(i+2)-un1(i+1)-un1(i)+un1(i-1))/(4d0*Deltax2)
              elseif (order ==  4) then
                 if (CIoption == 1) then
                    xd = 0.5d0 - alpham / Deltax ! same wether alpham is + or - 
                    fw = un1(i)
                    fe = un1(i+1)
                    fxw= fx(un1(i+1), un1(i-1), 2d0)
                    fxe= fx(un1(i+2), un1(i), 2d0)
                 elseif (CIoption == 2) then
                    xd = 1.5d0 - alpham / Deltax ! same wether alpham is + or - 
                    fw = un1(i-1)
                    fe = un1(i+2)
                    deno=2d0/3d0
                    fxw= fx(un1(i), un1(i-2), deno)
                    fxe= fx(un1(i+3), un1(i+1), deno)
                 endif
                 uinterp=cubic_interp (fw, fe, fxw,  fxe, xd)
              endif
              alpham=Deltat*uinterp
           enddo

!------------------------------------------------------------------------
! find hbef and Abef (initial position of particle at time level n-2)
!------------------------------------------------------------------------
        
           if (order .eq. 1) then
           
              hbef = hn2in(i) - ( hn2in(i+1) - hn2in(i-1) )*alpham / Deltax
              Abef = An2in(i) - ( An2in(i+1) - An2in(i-1) )*alpham / Deltax
              
           elseif (order .eq. 2) then
              
              hbef=hn2in(i) - (hn2in(i+1) - hn2in(i-1))*alpham / Deltax + &
                   2d0*(alpham**2d0)*(hn2in(i-1)-2d0*hn2in(i)+hn2in(i+1))/Deltax2

              Abef=An2in(i) - (An2in(i+1) - An2in(i-1))*alpham / Deltax + &
                   2d0*(alpham**2d0)*(An2in(i-1)-2d0*An2in(i)+An2in(i+1))/Deltax2
           
           elseif (order .eq. 4) then ! cubic interpolation

              if (alpham .ge. 0d0) then
                 xdn2= 1d0  - 2d0*alpham /Deltax
                 xdn1= 1d0  - alpham /Deltax
                 iw=i-1
                 ie=i
              else ! alpham .lt. 0d0
                 xdn2=-2d0*alpham/Deltax
                 xdn1=-1d0*alpham/Deltax
                 iw=i
                 ie=i+1
              endif
              
              fw=hn2in(iw)
              fe=hn2in(ie)
              fxw=fx(hn2in(iw+1), hn2in(iw-1), 2d0)
              fxe=fx(hn2in(ie+1), hn2in(ie-1), 2d0)
              hbef=cubic_interp (fw, fe, fxw, fxe, xdn2)
              fw=An2in(iw)
              fe=An2in(ie)
              fxw=fx(An2in(iw+1), An2in(iw-1), 2d0)
              fxe=fx(An2in(ie+1), An2in(ie-1), 2d0)
              Abef=cubic_interp (fw, fe, fxw, fxe, xdn2)

           endif

           if (SLlimiter) then ! LIMITER COULD BE IMPROVED FOR O3
             
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

           hbef = max(hbef, 0d0)
           Abef = max(Abef, 0d0)
           Abef = min(Abef, 1d0)

!------------------------------------------------------------------------                                        
! find right hand side terms rhsh and rhsA (time level n-1 = n1)
! minus sign added later in final calc of hout, Aout
!------------------------------------------------------------------------ 
          
           if (order .le. 2) then

              rhshprime=( hn1in(i+1)*(un1(i+2) - un1(i+1)) - &
                   hn1in(i-1)*(un1(i)   - un1(i-1)) ) / (2d0*Deltax2)
           
              rhsAprime=( An1in(i+1)*(un1(i+2) - un1(i+1)) - &
                   An1in(i-1)*(un1(i)   - un1(i-1)) ) / (2d0*Deltax2)

              if (order .eq. 1) then
                 rhsh = hn1in(i)*(un1(i+1)-un1(i))/Deltax - alpham*rhshprime
                 rhsA = An1in(i)*(un1(i+1)-un1(i))/Deltax - alpham*rhsAprime
                 
              elseif (order .eq. 2) then ! O2...O4 not coded yet
           
                 rhsh=hn1in(i)*(un1(i+1)-un1(i))/Deltax - alpham*rhshprime + &
                      (alpham**2d0)*(hn1in(i+1)*(un1(i+2)-un1(i+1)) - &
                      2d0*hn1in(i)*(un1(i+1)-un1(i)) + &
                      hn1in(i-1)*(un1(i)-un1(i-1)) ) / (2d0*(Deltax**3))
              
                 rhsA=An1in(i)*(un1(i+1)-un1(i))/Deltax - alpham*rhsAprime + &
                      (alpham**2d0)*(An1in(i+1)*(un1(i+2)-un1(i+1)) - &
                      2d0*An1in(i)*(un1(i+1)-un1(i)) + &
                      An1in(i-1)*(un1(i)-un1(i-1)) ) / (2d0*(Deltax**3))
              endif

           else
              fw=hn1in(iw)*div(iw)
              fe=hn1in(ie)*div(ie)
              fxw=fx(fe, hn1in(iw-1)*div(iw-1), 2d0)
              fxe=fx(hn1in(ie+1)*div(ie+1), fw, 2d0)
              rhsh=cubic_interp (fw, fe, fxw, fxe, xdn1)
              fw=An1in(iw)*div(iw)
              fe=An1in(ie)*div(ie)
              fxw=fx(fe, An1in(iw-1)*div(iw-1), 2d0)
              fxe=fx(An1in(ie+1)*div(ie+1), fw, 2d0)
              rhsA=cubic_interp (fw, fe, fxw, fxe, xdn1)
           endif

!------------------------------------------------------------------------
! find hout, Aout (after, time level n)
!------------------------------------------------------------------------ 

           hout(i) = hbef - 2d0*Deltat*rhsh
           Aout(i) = Abef - 2d0*Deltat*rhsA

        elseif (caseSL == 2) then ! upwind used close to walls

           flux=calc_flux(utp(i),utp(i+1),hn1in(i-1),hn1in(i), hn1in(i+1)) ! for h
           hout(i) = hn1in(i) - DtoverDx*flux
           flux=calc_flux(utp(i),utp(i+1),An1in(i-1),An1in(i), An1in(i+1)) ! for A                                                    
           Aout(i) = An1in(i) - DtoverDx*flux

        endif

        hout(i) = max(hout(i), 0d0)
        Aout(i) = max(Aout(i), 0d0)
        Aout(i) = min(Aout(i), 1d0)
        
     enddo

  endif
  
  return
end subroutine advection

subroutine calc_div (utp, div)
  use size
  use resolution  
  
  implicit none
  
  integer i
  double precision, intent(in) :: utp(1:nx+1)
  double precision, intent(out):: div(nx)

  do i = 1, nx
     div(i) = ( utp(i+1)-utp(i) ) / Deltax 
  enddo

end subroutine calc_div
  
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

function cubic_interp (fw, fe, fxw, fxe, xdist) result(finterp)
  use resolution

! see https://www.paulinternet.nl/?page=bicubic
! dist is between 0 (west=w) and 1 (east=e)

  double precision, intent(in) :: fw, fe, fxw, fxe, xdist    ! input
  double precision             :: a, b ! parameters for cubic interpolation
  double precision             :: finterp                    ! output
                                                                      
  a = 2d0*fw - 2d0*fe + fxw + fxe
  b = -3d0*fw + 3d0*fe - 2d0*fxw -fxe
  
  finterp = a*(xdist**3d0) + b*(xdist**2d0) + fxw*xdist + fw

end function cubic_interp

function calc_flux (ui, uip1, Tim1, Ti, Tip1) result(flux)
  use resolution

  double precision, intent(in) :: ui, uip1, Tim1, Ti, Tip1 ! T=tracer
  double precision :: flux, flux1, flux2

  if (ui .gt. Deltax/Deltat) print *, 'WARNING: u > dx/dt', ui

  if (ui .ge. 0d0) then ! left side of cell                                                                                     
     flux1 = ui*Tim1
  else
     flux1 = ui*Ti
  endif

  if (uip1 .ge. 0d0) then ! right side of cell                                                                                  
     flux2 = uip1*Ti
  else
     flux2 = uip1*Tip1
  endif

  flux = flux2 - flux1

end function calc_flux

function fx(fright, fleft, deno) result(dfdx)
      
      double precision, intent(in) :: fright, fleft                                 
      double precision, intent(in) :: deno
      double precision             :: dfdx ! output df/dx                                          

      dfdx = ( fright - fleft ) / deno
      
end function fx




