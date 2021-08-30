subroutine output_results(ts, expnb, solver, utp, zeta, eta)
  use size
  use resolution
  use global_var
  use shallow_water
  use MOMeqSW_output
  use rheology
  use option
  implicit none

  character filename*80

  integer :: i, k, Dt, Dx, adv
  integer, intent(in) :: ts, expnb, solver
  double precision, intent(in):: zeta(0:nx+1),eta(0:nx+1)
  double precision, intent(in)  :: utp(1:nx+1)
  double precision :: div(0:nx+1), sigma(0:nx+1), sig_norm(0:nx+1), zeta_norm(0:nx+1)
  double precision :: Erate(1:nx+1) ! KE rate loss/gain by rheology term
  character(Len=10) :: fldr

  if (adv_scheme .eq. 'upwind') then
    adv = 1
  elseif (adv_scheme .eq. 'upwindRK2') then
    adv = 2
  elseif (adv_scheme .eq. 'semilag') then
    adv = 3
  endif

  fldr='output.04/' !output folder

  Dt=int(Deltat) ! in s
  Dx=int(Deltax/1000d0) ! in km

  div(0) = 0d0
  div(nx+1) = 0d0
!  zeta(0) = 0d0
!  zeta(nx+1) = 0d0
  sigma(0) = 0d0
  sigma(nx+1) = 0d0
  sig_norm(0) = 0d0
  sig_norm(nx+1) = 0d0
  Erate = 0d0

  do i = 1, nx
     div(i) = (utp(i+1)-utp(i)) / Deltax ! calc divergence
     sigma(i) = (zeta(i)+eta(i))*div(i) - P_half(i)
     sig_norm(i) = (zeta(i)+eta(i))*div(i)*0.5d0/Pp_half(i) - 0.5d0*P_half(i)/Pp_half(i) ! norm by ice strength
!     zeta_norm(i) = zeta(i) / (zmax_par*Pp_half(i))
  enddo

  do i = 2, nx
    Erate(i) = utp(i) * ( sigma(i) - sigma(i-1) ) / Deltax
  enddo

  write (filename, '(A,"h_",i5.5,"s_",i3.3,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') fldr, Dt, &
        Dx,solver,IMEX, adv,BDF2,ts,expnb
  open (10, file = filename, status = 'unknown')

  write (filename, '(A,"A_",i5.5,"s_",i3.3,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') fldr, Dt, &
        Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (11, file = filename, status = 'unknown')

  write (filename, '(A,"u_",i5.5,"s_",i3.3,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') fldr, Dt, &
        Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (12, file = filename, status = 'unknown')

  write (filename, '(A,"div_",i5.5,"s_",i3.3,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') fldr, Dt, &
        Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (13, file = filename, status = 'unknown')

!  write (filename, '(A,"zeta_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_ts",i4.4,".",i2.2)') fldr, Dt,Dx, &
!       IMEX, adv,ts,expnb
!  open (14, file = filename, status = 'unknown')

!  write (filename, '(A,"sigma_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_ts",i4.4,".",i2.2)') fldr, Dt,Dx, &
!       IMEX, adv,ts,expnb
!  open (15, file = filename, status = 'unknown')

!  write (filename, '(A,"zeta_norm_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_ts",i4.4,".",i2.2)') fldr, Dt,Dx, &
!       IMEX, adv,ts,expnb
!  open (16, file = filename, status = 'unknown')

!  write (filename, '(A,"sig_norm_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_ts",i4.4,".",i2.2)') fldr, Dt,Dx, &
!       IMEX, adv,ts,expnb
!  open (17, file = filename, status = 'unknown')

  write (filename, '(A,"Er_",i5.5,"s_",i3.3,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') fldr, Dt, &
        Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (18, file = filename, status = 'unknown')


  write(10,10) ( h(i),       i = 0, nx+1 )
  write(11,10) ( A(i),       i = 0, nx+1 )
  write(13,10) ( div(i),     i = 0, nx+1 )
!  write(14,10) ( zeta(i),    i = 0, nx+1 )
!  write(15,10) ( sigma(i),   i = 0, nx+1 )
!  write(16,10) ( zeta_norm(i),    i = 0, nx+1 )
!  write(17,10) ( sig_norm(i),   i = 0, nx+1 )
  write(12,10) ( utp(i),       i = 1, nx+1 )
  write(18,10) ( Erate(i),     i = 1, nx+1 )

  do k = 10, 18
     close(k)
  enddo

  if (oceanSIM) then

  write (filename, '("output/etaw_",i5.5,"s_",i3.3,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') Dt, &
        Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (18, file = filename, status = 'unknown')

  write (filename, '("output/uw_",i5.5,"s_",i3.3,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') Dt, &
        Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (19, file = filename, status = 'unknown')

  write (filename, '("output/duwdt_",i5.5,"s_",i3.3,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') &
        Dt, Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (20, file = filename, status = 'unknown')

  write (filename, '("output/gdetawdx_",i5.5,"s_",i3.3,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') &
        Dt, Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (21, file = filename, status = 'unknown')

  write (filename, '("output/tauaw_",i5.5,"s_",i3.3,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') &
        Dt, Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (22, file = filename, status = 'unknown')

  write (filename, '("output/tauiw_",i5.5,"s_",i3.3,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') &
        Dt, Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (23, file = filename, status = 'unknown')

  write (filename, '("output/buw_",i5.5,"s_",i3.3,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') &
        Dt, Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (24, file = filename, status = 'unknown')


  write(18,10) ( etaw(i),      i = 0, nx+1 )
  write(19,10) ( uw(i),        i = 1, nx+1 )
  write(20,10) ( duwdt(i),     i = 1, nx+1 )
  write(21,10) ( gedetawdx(i), i = 1, nx+1 )
  write(22,10) ( tauaw(i),     i = 1, nx+1 )
  write(23,10) ( tauiw(i),     i = 1, nx+1 )
  write(24,10) ( buw(i),       i = 1, nx+1 )

  do k = 18, 24
     close(k)
  enddo

  endif

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

  write (filename, '("output/res_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_BDF2_",i1.1,"_ts",i4.4,"_k",i3.3,".",i2.2)') Dt, &
        Dx,IMEX,adv,BDF2,ts,k,expnb
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

  write (filename, '("output/Nbite_",i5.5,"s_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_BDF2_",i1.1,".",i2.2)') Dt, &
   Dx,IMEX,adv,BDF2,expnb
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

  write (filename, '("output/iniL2norm_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_BDF2_",i1.1,".",i2.2)') Dt,&
   Dx,IMEX,adv,BDF2,expnb
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

  write (filename, '("output/uk1_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_BDF2_",i1.1,"_ts",i4.4,"_k",i3.3,".dat")') Dt,&
        Dx,IMEX, adv,BDF2,ts,k
  open (10, file = filename, status = 'unknown')

  write (filename, '("output/du_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_BDF2_",i1.1"_ts",i4.4,"_k",i3.3,".dat")') Dt,&
        Dx,IMEX, adv,BDF2,ts,k
  open (11, file = filename, status = 'unknown')

  write(10,10) ( utp(i),       i = 1, nx+1 )
  write(11,10) ( du(i),       i = 1, nx+1 )

  close(10)
  close(11)

10 format (1x, 1000(f25.20, 1x))

  return
end subroutine output_u_and_du

subroutine output_diag_stress(ts, expnb, idiag)

! output diagnostic (ice-ocean and vice versa) stress at i=idiag

  use resolution
  use diag_stress
  use option
  implicit none

  character filename*70

  integer :: Dt, Dx, adv
  integer, intent(in) :: ts, expnb, idiag
  double precision :: ratio

  if (adv_scheme .eq. 'upwind') then
    adv = 1
  elseif (adv_scheme .eq. 'upwindRK2') then
    adv = 2
  endif

  Dt=int(Deltat/60d0) ! in min
  Dx=int(Deltax/1000d0) ! in km

  ratio = tauwidiag/tauiwdiag

  write (filename, '("output/diag_stress_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_BDF2_",i1.1,".",i2.2)') Dt,&
   Dx,IMEX,adv,BDF2,expnb
  open (10, file = filename, access = 'append')

  write(10,10) ts, idiag, tauwidiag, tauiwdiag, ratio, tauaidiag

  close(10)

10 format (i5,1x, i5,1x,f12.8,1x,f12.8,1x,f12.8,1x,f12.8)

  return
end subroutine output_diag_stress
