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
  double precision :: zeta(0:nx+1), eta(0:nx+1)
  double precision :: Cw(1:nx+1), Fpos(1:nx+1)
  double precision :: upos(1:nx+1), b(1:nx+1)

!  double precision xpos(nvar), xneg(nvar), x(nvar),rhs(nvar)

!  double precision Fpos(nvar),Fneg(nvar)
!  double precision epsilon,v(nvar),Jv(nvar)
      
  do i=1, nx+1

     upos(i) = uk1(i) + epsilon * v(i)
!     xneg(i) = x(i) - epsilon * v(i)

  enddo
        
  if (IMEX .eq. 2) then ! IMEX method 2 only (WATCHOUT hpos for precond...)
    call advection (upts, upos, hpts, Apts, h, A) ! advection scheme for tracers
    call ice_strength () ! Pp_half is Pp/2 where Pp is the ice strength
  endif
  
  call viscouscoefficient (upos, zeta, eta)
  call bvect (tauair, upts, b)
  call Cw_coefficient (upos, Cw)
  call Fu (upos, zeta, eta, Cw, b, Fpos)

  do i=1, nx+1

     Jv(i) = ( Fpos(i)-F_uk1(i) ) / epsilon
            
  enddo

  return
end subroutine JacfreeVec
      
      
subroutine formJacobian (utp, Futp, upts, tauair, ts, k)

  use size
  use global_var
  use option
  
  implicit none
  
  integer :: i
  integer, intent(in) :: ts, k
  double precision, intent(in) :: Futp(1:nx+1), utp(1:nx+1)
  double precision, intent(in) :: upts(1:nx+1), tauair(1:nx+1)
  double precision :: zeta(0:nx+1), eta(0:nx+1), epsilon
  double precision :: Cw(1:nx+1), Fpos(1:nx+1)
  double precision :: uele(1:nx+1), upos(1:nx+1), b(1:nx+1)
  
  double precision :: Jleft(1:nx+1), J(1:nx+1), Jright(1:nx+1)

  epsilon=1d-10
  
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
  call bvect (tauair, upts, b)
  call Cw_coefficient (upos, Cw)
  call Fu (upos, zeta, eta, Cw, b, Fpos)
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
  call bvect (tauair, upts, b)
  call Cw_coefficient (upos, Cw)
  call Fu (upos, zeta, eta, Cw, b, Fpos)
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
  call bvect (tauair, upts, b)
  call Cw_coefficient (upos, Cw)
  call Fu (upos, zeta, eta, Cw, b, Fpos)
  Jright(i)=(Fpos(i)-Futp(i))/epsilon
  else
  Jright(i)=0d0
  endif

  uele(i+1)=0d0

  print *, 'Jacobian', ts, k, i, Jleft(i), J(i), Jright(i)

  enddo

  return
end subroutine formJacobian

subroutine formA (utp, zeta, eta, Cw, ts, k)
  use size
  use resolution
  use properties
  use global_var
  use option

  implicit none
  
  integer, intent(in) :: ts, k
  integer :: i

  double precision, intent(in)  :: utp(1:nx+1)
  double precision, intent(in)  :: zeta(0:nx+1), eta(0:nx+1)
  double precision, intent(in)  :: Cw(1:nx+1)
  double precision :: Aleft(1:nx+1), Adiag(1:nx+1), Aright(1:nx+1)

  double precision :: h_at_u

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

  enddo

  return
end subroutine formA
      


