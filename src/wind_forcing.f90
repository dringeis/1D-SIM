! defined at the u location

subroutine wind_forcing (tauair, ts)

  use size
  use forcing
  use resolution
  use option

  implicit none

  integer :: i, nf
  integer, intent(in) :: ts
  double precision :: speed, period, modulation, pi, apar, tau_f
  double precision, intent(out) :: tauair(1:nx+1) ! air drag

  nf = 129
  tau_f = 0.5
  speed = 10d0 ! [m/s]
  period = 6d0*3600d0 ! period of cos in seconds (set to 3 days)
  pi = 3.14159265d0
  apar = 6d0*3600d0

  tauair(1)    = 0d0 ! close bc
  tauair(nx+1) = 0d0 ! close bc

    if (constant_wind) then

      do i = 2, nf ! with apar = 6*3600, tauair is (1-e^-2) after 12 hours.
        ! tauair(i) = Cda * (speed)**2d0
        !	tauair(i) = (Cda * (speed)**2d0)*(1d0-exp(-1d0*ts*Deltat/apar)) ! at n
        ! print *, 'tauair', ts, 100d0*tauair(50)/(Cda * (speed)**2d0)
        tauair(i) = 0
      enddo
      do i = nf+1, nx
        tauair(i) = tau_f
      enddo

    else
      modulation = sin(2*pi*ts*Deltat/period)

      do i = 2, nx
        tauair(i) = modulation * Cda * (speed)**2d0
      enddo

    endif

  return
end subroutine wind_forcing






