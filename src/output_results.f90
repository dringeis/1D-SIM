subroutine output_results(ts, expnb, zeta, eta)
  use size
  use resolution
  use global_var
  use rheology
  implicit none

  character filename*30

  integer :: i, k
  integer, intent(in) :: ts, expnb
  double precision, intent(in):: zeta(0:nx+1),eta(0:nx+1)
  double precision :: div(0:nx+1), sigma(0:nx+1), sig_norm(0:nx+1), zeta_norm(0:nx+1)

  div(0) = 0d0
  div(nx+1) = 0d0

  do i = 1, nx
     div(i) = (u(i+1)-u(i)) / Deltax ! calc divergence
     sigma(i) = (zeta(i)+eta(i))*div(i) - P_half(i)
     sig_norm(i) = (zeta(i)+eta(i))*div(i)*0.5d0/P_half(i) - 0.5d0
     zeta_norm(i) = zeta(i) / (zmax_par*Pp_half(i))
  enddo
  
  write (filename, '("output/h_",i3.3,".",i2.2)') ts,expnb
  open (10, file = filename, status = 'unknown')
  
  write (filename, '("output/A_",i3.3,".",i2.2)') ts,expnb
  open (11, file = filename, status = 'unknown')

  write (filename, '("output/u_",i3.3,".",i2.2)') ts,expnb
  open (12, file = filename, status = 'unknown')

  write (filename, '("output/div_",i3.3,".",i2.2)') ts,expnb
  open (13, file = filename, status = 'unknown')

  write (filename, '("output/zeta_",i3.3,".",i2.2)') ts,expnb
  open (14, file = filename, status = 'unknown')

  write (filename, '("output/sigma_",i3.3,".",i2.2)') ts,expnb
  open (15, file = filename, status = 'unknown')

  write (filename, '("output/zeta_norm",i3.3,".",i2.2)') ts,expnb
  open (16, file = filename, status = 'unknown')

  write (filename, '("output/sig_norm",i3.3,".",i2.2)') ts,expnb
  open (17, file = filename, status = 'unknown')


  write(10,10) ( h(i),       i = 0, nx+1 )
  write(11,10) ( A(i),       i = 0, nx+1 )
  write(13,10) ( div(i),     i = 0, nx+1 )
  write(14,10) ( zeta(i),    i = 0, nx+1 )
  write(15,10) ( sigma(i),   i = 0, nx+1 )
  write(16,10) ( zeta_norm(i),    i = 0, nx+1 )
  write(17,10) ( sig_norm(i),   i = 0, nx+1 )
  write(12,10) ( u(i),       i = 1, nx+1 )

  do k = 10, 17
     close(k)
  enddo

10 format (1x, 1000(f25.20, 1x))

  return
end subroutine output_results

subroutine output_residual(ts, k, expnb, F)
  use size

  implicit none

  character filename*40

  integer :: i
  integer, intent(in) :: ts, k, expnb
  double precision, intent(in):: F(1:nx+1)

  write (filename, '("output/residual_",i4.4,"_",i4.4,".",i2.2)') ts,k,expnb
  open (11, file = filename, status = 'unknown')

  write(11,10) ( F(i), i = 1, nx+1 )

  close(11)
  
10 format (1x, 1000(f25.18, 1x))

  return
end subroutine output_residual

subroutine output_nb_Newton_or_OL_ite(k)
  use resolution
  use option
  implicit none

  character filename*35

  integer :: Dt
  integer, intent(in) :: k

  Dt=int(Deltat)

  write (filename, '("output/Newton_ite_",i5.5".",i2.2)') Dt,IMEX
  open (10, file = filename, access = 'append')
  
  write(10,10) k-1

  close(10)

10 format (i4)

  return
end subroutine output_nb_Newton_or_OL_ite
