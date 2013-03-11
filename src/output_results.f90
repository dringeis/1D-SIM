subroutine output_results(ts, expnb, utp, zeta, eta)
  use size
  use resolution
  use global_var
  use rheology
  use option
  implicit none

  character filename*60

  integer :: i, k, Dt, Dx, adv
  integer, intent(in) :: ts, expnb
  double precision, intent(in):: zeta(0:nx+1),eta(0:nx+1)
  double precision, intent(in)  :: utp(1:nx+1)
  double precision :: div(0:nx+1), sigma(0:nx+1), sig_norm(0:nx+1), zeta_norm(0:nx+1)

  if (adv_scheme .eq. 'upwind') then
    adv = 1
  elseif (adv_scheme .eq. 'upwindRK2') then
    adv = 2
  endif

  Dt=int(Deltat) ! in s
  Dx=int(Deltax/1000d0) ! in km

  div(0) = 0d0
  div(nx+1) = 0d0

  do i = 1, nx
     div(i) = (utp(i+1)-utp(i)) / Deltax ! calc divergence
     sigma(i) = (zeta(i)+eta(i))*div(i) - P_half(i)
     sig_norm(i) = (zeta(i)+eta(i))*div(i)*0.5d0/P_half(i) - 0.5d0
     zeta_norm(i) = zeta(i) / (zmax_par*Pp_half(i))
  enddo
  
  write (filename, '("output/h_",i5.5,"s_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_CN",i1.1,"_ts",i6.6,".",i2.2)') Dt,Dx, &
		    IMEX, adv,CN,ts,expnb
  open (10, file = filename, status = 'unknown')
  
  write (filename, '("output/A_",i5.5,"s_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_CN",i1.1,"_ts",i6.6,".",i2.2)') Dt,Dx, &
		    IMEX, adv,CN,ts,expnb
  open (11, file = filename, status = 'unknown')

  write (filename, '("output/u_",i5.5,"s_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_CN",i1.1,"_ts",i6.6,".",i2.2)') Dt,Dx, &
		    IMEX, adv,CN,ts,expnb
  open (12, file = filename, status = 'unknown')

  write (filename, '("output/div_",i5.5,"s_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_CN",i1.1,"_ts",i6.6,".",i2.2)') Dt,Dx, &
		    IMEX, adv,CN,ts,expnb
  open (13, file = filename, status = 'unknown')

!  write (filename, '("output/zeta_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_ts",i4.4,".",i2.2)') Dt,Dx, &
!		    IMEX, adv,ts,expnb
!  open (14, file = filename, status = 'unknown')

!  write (filename, '("output/sigma_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_ts",i4.4,".",i2.2)') Dt,Dx, &
!		    IMEX, adv,ts,expnb
!  open (15, file = filename, status = 'unknown')

!  write (filename, '("output/zeta_norm_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_ts",i4.4,".",i2.2)') Dt,Dx, &
!		    IMEX, adv,ts,expnb
!  open (16, file = filename, status = 'unknown')

!  write (filename, '("output/sig_norm_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_ts",i4.4,".",i2.2)') Dt,Dx, &
!		    IMEX, adv,ts,expnb
!  open (17, file = filename, status = 'unknown')


  write(10,10) ( h(i),       i = 0, nx+1 )
  write(11,10) ( A(i),       i = 0, nx+1 )
  write(13,10) ( div(i),     i = 0, nx+1 )
!  write(14,10) ( zeta(i),    i = 0, nx+1 )
!  write(15,10) ( sigma(i),   i = 0, nx+1 )
!  write(16,10) ( zeta_norm(i),    i = 0, nx+1 )
!  write(17,10) ( sig_norm(i),   i = 0, nx+1 )
  write(12,10) ( utp(i),       i = 1, nx+1 )

  do k = 10, 17
     close(k)
  enddo

10 format (1x, 1000(f25.20, 1x))

  return
end subroutine output_results

subroutine output_residual(ts, k, expnb, F)
  use size
  use resolution
  use option

  implicit none

  character filename*60

  integer :: i, Dt, Dx, adv
  integer, intent(in) :: ts, k, expnb
  double precision, intent(in):: F(1:nx+1)
 
  if (adv_scheme .eq. 'upwind') then
    adv = 1
  elseif (adv_scheme .eq. 'upwindRK2') then
    adv = 2
  endif

  Dt=int(Deltat/60d0) ! in min
  Dx=int(Deltax/1000d0) ! in km
  
  write (filename, '("output/res_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_CN",i1.1,"_ts",i4.4,"_k",i3.3,".",i2.2)') Dt,Dx,IMEX, &
		    adv,CN,ts,k,expnb
  open (11, file = filename, status = 'unknown')
  
  write(11,10) ( F(i), i = 1, nx+1 )

  close(11)
  
10 format (1x, 1000(f25.18, 1x))

  return
end subroutine output_residual

subroutine output_nb_ite(ts, k, fgmres_per_ts, expnb)
  use resolution
  use option
  implicit none

  character filename*60

  integer, intent(in) :: ts, k, fgmres_per_ts, expnb
  integer :: Dt, Dx, adv

  if (adv_scheme .eq. 'upwind') then
    adv = 1
  elseif (adv_scheme .eq. 'upwindRK2') then
    adv = 2
  endif

  Dt=int(Deltat) ! in s
  Dx=int(Deltax/1000d0) ! in km

  write (filename, '("output/Nbite_",i5.5,"s_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_CN",i1.1,".",i2.2)') Dt,Dx,&
	 IMEX,adv,CN,expnb
  open (10, file = filename, access = 'append')
  
  write(10,10) ts,k-1,fgmres_per_ts

  close(10)

10 format (i5,1x,i4,1x,i5)

  return
end subroutine output_nb_ite

subroutine output_ini_L2norm(ts, L2norm, expnb)
  use resolution
  use option
  implicit none

  character filename*60

  integer :: Dt, Dx, adv
  integer, intent(in) :: ts, expnb
  double precision, intent(in) :: L2norm

  if (adv_scheme .eq. 'upwind') then
    adv = 1
  elseif (adv_scheme .eq. 'upwindRK2') then
    adv = 2
  endif

  Dt=int(Deltat/60d0) ! in min
  Dx=int(Deltax/1000d0) ! in km

  write (filename, '("output/iniL2norm_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_CN",i1.1,".",i2.2)') Dt,Dx,&
	 IMEX,adv,CN,expnb
  open (10, file = filename, access = 'append')
  
  write(10,10) ts,L2norm

  close(10)

10 format (i5,1x,f15.12)

  return
end subroutine output_ini_L2norm

subroutine output_u_and_du ( ts, k, utp, du )
  use size
  use resolution
!  use global_var
  use rheology
  use option
  implicit none

  character filename*60

  integer :: i, Dt, Dx, adv
  integer, intent(in) :: ts, k
  double precision, intent(in)  :: utp(1:nx+1), du(1:nx+1)

  if (adv_scheme .eq. 'upwind') then
    adv = 1
  elseif (adv_scheme .eq. 'upwindRK2') then
    adv = 2
  endif

  Dt=int(Deltat/60d0) ! in min
  Dx=int(Deltax/1000d0) ! in km

  write (filename, '("output/uk1_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_CN",i1.1,"_ts",i4.4,"_k",i3.3,".dat")') Dt,Dx, &
		    IMEX, adv,CN,ts,k
  open (10, file = filename, status = 'unknown')
  
  write (filename, '("output/du_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,,"_CN",i1.1"_ts",i4.4,"_k",i3.3,".dat")') Dt,Dx, &
		    IMEX, adv,CN,ts,k
  open (11, file = filename, status = 'unknown')
  
  write(10,10) ( utp(i),       i = 1, nx+1 )
  write(11,10) ( du(i),       i = 1, nx+1 )

  close(10)
  close(11)

10 format (1x, 1000(f25.20, 1x))

  return
end subroutine output_u_and_du