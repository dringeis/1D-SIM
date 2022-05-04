!****************************************************************************
!     initial thickness and velocity field (u is at t=0, h is at t=Deltat/2)
!****************************************************************************

! no restart possible for the moment jfl WATCHOUT
! no open bcs possible for the moment jfl WATCHOUT

subroutine ini_get (utp, restart, expres, ts_res)

  use size
  use global_var
  use shallow_water
  use MOMeqSW_output
  use option

  implicit none

  logical, intent(in) :: restart
  integer, intent(in) :: expres, ts_res
  integer :: i
  double precision, intent(inout)  :: utp(1:nx+1)
  double precision :: rdnb, small

  character(LEN=30) filename  ! restart file name

  small = 0.0001d0

  allocate(etaw(0:nx+1), etawn1(0:nx+1), etawn2(0:nx+1))
  allocate(uw(1:nx+1), uwn1(1:nx+1), uwn2(1:nx+1))

  if (oceanSIM) then
     allocate(duwdt(1:nx+1), gedetawdx(1:nx+1), buw(1:nx+1))
     allocate(tauiw(1:nx+1), tauaw(1:nx+1))
  endif

  !utp(1)    = 0d0 ! close bc
  !utp(nx+1) = 0d0 ! close bc
  !h(0)    = 0d0
  !h(nx+1) = 0d0
  !A(0)    = 0d0
  !A(nx+1) = 0d0
  uw      = 0d0
  uwn1    = 0d0
  uwn2    = 0d0
  etaw    = 0d0
  etawn1  = 0d0
  etawn2  = 0d0

  scaling=1d0 ! initialize scaling field (only used for JFNK)

  if (restart) then
     print *, 'Restart code should be verified'
     stop
     write (filename,'("output/h_",i3.3,".",i2.2)') ts_res, expres
     open (10, file = filename, status = 'old')

     write (filename,'("output/A_",i3.3,".",i2.2)') ts_res, expres
     open (11, file = filename, status = 'old')

     write (filename,'("output/u_",i3.3,".",i2.2)') ts_res, expres
     open (12, file = filename, status = 'old')

     read (10,*) ( h(i), i = 0, nx+1 )
     read (11,*) ( A(i), i = 0, nx+1 )
     read (12,*) ( utp(i), i = 1, nx+1 )

     close(10)
     close(11)
     close(12)

  else ! specify initial fields

  do i = 1, nx+1
!     call random_number(rdnb)
!     u(i) = small*(rdnb-0.5d0) !small random nb added to 1st initial guess
     utp(i) = 0d0
     uw(i)  = 0d0
  enddo

  do i = 0, nx+1
     ! h(i) = 0d0
     ! A(i) = 0d0
!     A(i) = i/(nx*1d0) - 0.5d0/(1d0*nx) ! 0 at West wall and 1 at East wall
!     h(i) = max(1d-06, h(i))
     bathy(i)=200d0
  enddo

  do i = 0, nx+1
     h(i) = 1d0
     A(i) = 1d0
!     A(i) = i/(nx*1d0) - 0.5d0/(1d0*nx) ! 0 at West wall and 1 at East wall
!     h(i) = max(1d-06, h(i))
!     bathy(i)=100d0
  enddo

  if (oceanSIM) then
     uwn1=uw
     uwn2=uw
     etawn1=etaw
     etawn2=etaw
     duwdt=0d0
     gedetawdx=0d0
     tauiw=0d0
     tauaw=0d0
     buw=0d0
  endif

!  do i = 11, nx-10
!     h(i) = 1d0
!     A(i) = 0.7d0
!     A(i) = i/(nx*1d0) - 0.5d0/(1d0*nx) ! 0 at West wall and 1 at East wall
!     h(i) = max(1d-06, h(i))
!     bathy(i)=100d0
!  enddo

!  do i = 1, 20
!     bathy(i)=10d0
!  enddo

!  do i = 1, 19
!     h(i) = 1d0
!     A(i) = 1d0
!  enddo

  endif

  return
end subroutine ini_get






