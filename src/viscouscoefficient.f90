subroutine viscouscoefficient(utp, zeta, eta)
  use size
  use resolution
  use rheology
  use option
  use global_var

  implicit none

  integer :: i

  double precision, intent(in) :: utp(1:nx+1)
  double precision, intent(out):: zeta(0:nx+1), eta(0:nx+1)
  double precision :: dudx, deno, denonum, denomin

  denomin=2d-09

  do i = 1, nx ! for tracer points

     dudx = ( utp(i+1) - utp(i) ) / Deltax

    if ( regularization .eq. 'tanh' ) then

     deno = alpha*sqrt( (dudx)**2d0 + small2 ) ! small2 is there to avoid div by zero
!     deno = alpha*(abs(dudx))
!     deno = max( deno, 1d-30 )
     zeta(i) = ((Pp_half(i)+Tp_half(i))/denomin)*tanh(denomin*(1d0/deno))

    elseif ( regularization .eq. 'Kreyscher' ) then

      deno = alpha*sqrt( (dudx)**2d0)
      zeta(i) = (Pp_half(i)+Tp_half(i)) / ( deno + denomin )

    elseif ( regularization .eq. 'capping' ) then

      deno = max((alpha*sqrt( (dudx)**2d0)), denomin)
      zeta(i) = (Pp_half(i)+Tp_half(i)) / deno

    else
      print *, 'Wrong regularization'
      stop

    endif

    eta(i)  = zeta(i) * e_2

    if (rep_closure) then  ! replacement closure (Kreysher et al. 2000)

      if ( regularization .eq. 'tanh' ) then
        P_half(i) = (Pp_half(i)-Tp_half(i)) * ( deno / denomin ) * tanh(denomin*(1d0/deno))
      elseif ( regularization .eq. 'Kreyscher' ) then
        P_half(i) = (Pp_half(i)-Tp_half(i)) * deno / ( deno + denomin )
      elseif ( regularization .eq. 'capping' ) then
        denonum=alpha*sqrt( (dudx)**2d0)
        P_half(i) = (Pp_half(i)-Tp_half(i)) * denonum / deno
      endif

    else

      P_half(i) = Pp_half(i)-Tp_half(i) ! P_half includes tensile strength

    endif

  enddo

  if (linear_viscous) then
     zeta = 1d08
     eta  = e_2 * zeta
     P_half = 0d0
  endif

  return
end subroutine viscousCoefficient



