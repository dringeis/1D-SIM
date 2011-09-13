!****************************************************************************
!     initial thickness and velocity field (u is at t=0, h is at t=Deltat/2)
!****************************************************************************

! no restart possible for the moment jfl WATCHOUT
! no open bcs possible for the moment jfl WATCHOUT

subroutine ini_get (upts)

  use size
  use global_var

  implicit none
     
  integer :: i
  double precision :: rdnb, small
  double precision, intent(out) :: upts(1:nx+1) ! u previous time step

  small = 0.0001d0

  upts(:) = 0d0 ! because no restart

  u(1)    = 0d0 ! close bc
  u(nx+1) = 0d0 ! close bc
  h(0)    = 0d0
  h(nx+1) = 0d0
  A(0)    = 0d0
  A(nx+1) = 0d0

  do i = 2, nx
!     call random_number(rdnb)
!     u(i) = small*(rdnb-0.5d0) !small random nb added to 1st initial guess  
     u(i) = 0d0
  enddo

  do i = 1, nx
     h(i) = 0.5d0
     A(i) = 0.95d0
!     A(i) = i/(nx*1d0) - 0.5d0/(1d0*nx) ! 0 at West wall and 1 at East wall
  enddo

  return
end subroutine ini_get






