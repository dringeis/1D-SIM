subroutine JacfreeVec (v, Jv, F_uk1, uk1, upts, tauair, epsilon)

  use size
  use global_var
  use option
  
  implicit none
  
  integer :: i
  double precision, intent(in) :: v(1:nx+1), F_uk1(1:nx+1), uk1(1:nx+1)
  double precision, intent(in) :: upts(1:nx+1), tauair(1:nx+1)
  double precision, intent(out):: Jv(1:nx+1)
  double precision, intent(in) :: epsilon
  double precision :: zeta(0:nx+1), eta(0:nx+1), hmidp(0:nx+1), Amidp(0:nx+1) ! p = pos
  double precision :: umidp(1:nx+1), Cw(1:nx+1), Fpos(1:nx+1), Rpos(1:nx+1)
  double precision :: upos(1:nx+1)!, b(1:nx+1)

!  double precision xpos(nvar), xneg(nvar), x(nvar),rhs(nvar)

!  double precision Fpos(nvar),Fneg(nvar)
!  double precision epsilon,v(nvar),Jv(nvar)
      
  do i=1, nx+1

     upos(i) = uk1(i) + epsilon * v(i)
!     xneg(i) = x(i) - epsilon * v(i)

  enddo
        
  
  if (IMEX .eq. 2) then ! IMEX method 1 or 2
     call advection (upts, upos, hpts, Apts, h, A) ! advection scheme for tracers
     if ( CN .eq. 0 ) then
	call ice_strength (h, A) ! Pp_half is Pp/2 where Pp is the ice strength
     elseif ( CN .eq. 1 ) then
	hmidp=(h + hpts)/2d0 ! on pourait avoir h=(h+hpts)/2 ???
	Amidp=(A + Apts)/2d0
	call ice_strength (hmidp, Amidp)
     endif  
  endif
  
  if ( CN .eq. 0 ) then
	  call viscouscoefficient (upos, zeta, eta) ! u is u^k-1
	  call Cw_coefficient (upos, Cw)            ! u is u^k-1
	  call calc_R (upos, zeta, eta, Cw, tauair, Rpos)
	  call Fu (upos, upts, h, Rpos, Fpos) 
  elseif ( CN .eq. 1 ) then
	  umidp=(upos + upts)/2d0
	  call viscouscoefficient (umidp, zeta, eta) ! u is u^k-1
	  call Cw_coefficient (umidp, Cw)
	  call calc_R (umidp, zeta, eta, Cw, tauair, Rpos)
	  call Fu (upos, upts, hmidp, Rpos, Fpos)
  endif

  do i=1, nx+1

     Jv(i) = ( Fpos(i)-F_uk1(i) ) / epsilon
            
  enddo

  return
end subroutine JacfreeVec
      
      
subroutine formJacobian (utp, Futp, upts, tauair, ts, k, Jleft, J, Jright)

  use size
  use global_var
  use resolution
  use option
  
  implicit none
  
  character filename*60
  integer, intent(in) :: ts, k
  integer :: i, Dt, Dx, adv
  
  double precision, intent(in) :: Futp(1:nx+1), utp(1:nx+1)
  double precision, intent(in) :: upts(1:nx+1), tauair(1:nx+1)
  double precision :: zeta(0:nx+1), eta(0:nx+1), epsilon
  double precision :: Cw(1:nx+1), Fpos(1:nx+1)
  double precision :: uele(1:nx+1), upos(1:nx+1), b(1:nx+1)
  
  double precision, intent(out) :: Jleft(1:nx+1), J(1:nx+1), Jright(1:nx+1)

  Dt=int(Deltat/60d0) ! in min
  Dx=int(Deltax/1000d0) ! in km

  if (adv_scheme .eq. 'upwind') then
    adv = 1
  elseif (adv_scheme .eq. 'upwindRK2') then
    adv = 2
  endif

  epsilon=1d-10

  Jleft = 0d0
  J     = 0d0
  Jright= 0d0

  do i = 2, nx
  
  uele=0d0  
  
!------- left: i-1 --------    
  
  if (i .gt. 2) then
  uele(i-1)=epsilon
  upos = utp + uele
  
  if (IMEX .eq. 2) then ! IMEX method 2 only (WATCHOUT hpos for precond...)
    call advection (upts, upos, hpts, Apts, h, A) ! advection scheme for tracers
    call ice_strength () ! Pp_half is Pp/2 where Pp is the ice strength
  endif
  call viscouscoefficient (upos, zeta, eta)
  stop ! (need to modify or latest changes)
!  call bvect (tauair, upts, b)
  call Cw_coefficient (upos, Cw)
!  call Fu (upos, zeta, eta, Cw, b, Fpos)
  Jleft(i)=(Fpos(i)-Futp(i))/epsilon
  else
  Jleft(i)=0d0
  endif
  
  uele(i-1)=0d0
  
!------- diagonal: i --------    

  uele(i)=epsilon
  upos = utp + uele
  
  if (IMEX .eq. 2) then ! IMEX method 2 only (WATCHOUT hpos for precond...)
    call advection (upts, upos, hpts, Apts, h, A) ! advection scheme for tracers
    call ice_strength () ! Pp_half is Pp/2 where Pp is the ice strength
  endif
  call viscouscoefficient (upos, zeta, eta)
!  call bvect (tauair, upts, b)
  call Cw_coefficient (upos, Cw)
!  call Fu (upos, zeta, eta, Cw, b, Fpos)
  J(i)=(Fpos(i)-Futp(i))/epsilon
  
  uele(i)=0d0

!------- right: i+1 --------    
    
  if (i .lt. nx) then
  uele(i+1)=epsilon      
  upos = utp + uele
  
  if (IMEX .eq. 2) then ! IMEX method 2 only (WATCHOUT hpos for precond...)
    call advection (upts, upos, hpts, Apts, h, A) ! advection scheme for tracers
    call ice_strength () ! Pp_half is Pp/2 where Pp is the ice strength
  endif
  call viscouscoefficient (upos, zeta, eta)
!  call bvect (tauair, upts, b)
  call Cw_coefficient (upos, Cw)
!  call Fu (upos, zeta, eta, Cw, b, Fpos)
  Jright(i)=(Fpos(i)-Futp(i))/epsilon
  else
  Jright(i)=0d0
  endif

  uele(i+1)=0d0

!  print *, i, Jleft(i), J(i), Jright(i)

  enddo

!  write (filename, '("output/J_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_ts",i4.4,"_k",i3.3,".dat")') Dt,Dx,IMEX, &
!		    adv,ts,k
!  open (11, file = filename, status = 'unknown')
  
!  do i=1,nx+1
!    write(11,10) Jleft(i), J(i), Jright(i)
!  enddo
  
!  close(11)
  
!  10 format (1x,f25.15, 1x, f25.15, 1x, f25.15)

  return
end subroutine formJacobian

subroutine formA (utp, zeta, eta, Cw, ts, k, Aleft, Adiag, Aright)
  use size
  use resolution
  use properties
  use global_var
  use resolution
  use option

  implicit none
  
  character filename*60
  integer, intent(in) :: ts, k
  integer :: i, Dt, Dx, adv

  double precision, intent(in)  :: utp(1:nx+1)
  double precision, intent(in)  :: zeta(0:nx+1), eta(0:nx+1)
  double precision, intent(in)  :: Cw(1:nx+1)
  double precision, intent(out) :: Aleft(1:nx+1), Adiag(1:nx+1), Aright(1:nx+1)

  double precision :: h_at_u
 
  Dt=int(Deltat/60d0) ! in min
  Dx=int(Deltax/1000d0) ! in km
  
  if (adv_scheme .eq. 'upwind') then
    adv = 1
  elseif (adv_scheme .eq. 'upwindRK2') then
    adv = 2
  endif
  
  Aleft = 0d0
  Adiag = 0d0
  Aright = 0d0

  do i = 2, nx

!------------------------------------------------------------------------
!    rhoice*h*du/dt : tendency term, advection of momentum is neglected
!------------------------------------------------------------------------

     h_at_u = ( h(i) + h(i-1) ) / 2d0

     Adiag(i) = ( rho * h_at_u ) / Deltat
     
!------------------------------------------------------------------------
!     Cw*u : water drag term
!------------------------------------------------------------------------
     
     Adiag(i) = Adiag(i) + Cw(i)
     
!------------------------------------------------------------------------
!     -d ( (zeta+eta) du/dx ) / dx : rheology term
!------------------------------------------------------------------------

     Aleft(i) = - (zeta(i-1)+eta(i-1)) / Deltax2
     Adiag(i) = Adiag(i) + (zeta(i)+eta(i)+zeta(i-1)+eta(i-1)) / Deltax2
     Aright(i)= - (zeta(i)+eta(i)) / Deltax2

!     print *, i, Aleft(i), Adiag(i), Aright(i)

  enddo

!  write (filename, '("output/Amat_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_ts",i4.4,"_k",i3.3,".dat")') Dt,Dx,IMEX, &
!		    adv,ts,k
!  open (11, file = filename, status = 'unknown')
  
!  do i=1,nx+1
!    write(11,10) Aleft(i), Adiag(i), Aright(i)
!  enddo
  
!  close(11)
  
!  10 format (1x,f25.15, 1x, f25.15, 1x, f25.15)

  return
end subroutine formA
      


