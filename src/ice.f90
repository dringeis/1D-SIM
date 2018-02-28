!*************************************************************************
!     program ice:
!       1D model that calculate the ice thickness (h), concentration (A) 
!       and ice velocity (u).
!       
!       momentum equation:  rho*h(u^n-u^n-1)/Deltat = f(u^n,h^n-1,A^n-1)
!                           This equation is solved implicitly for u^n
!                           (u). u^n-1 is the previous time step solution.
!
!       continuity equation:h^n = f(h^n-1, u^n)  
!                           The new value of h (and A) is obtained by 
!                           advecting h^n-1 with u^n. 
!       
!       author: JF Lemieux
!       version: 1.0 (20 april 2012)
!
!************************************************************************

program ice

  use size
  use rheology
  use forcing
  use properties
  use resolution
  use global_var
  use shallow_water
  use MOMeqSW_output
  use numerical
  use option
  
  implicit none

  logical :: p_flag, restart
  integer :: i, ii, ts, tsini, nstep, tsfin, k, s, Nmax_OL, solver, idiag
  integer :: out_step(5), expnb, expres, ts_res, fgmres_its, fgmres_per_ts
  integer, save :: Nfail, meanN ! nb of failures, mean Newton ite per ts
  double precision :: e, rhoair, Cdair, Cdwater
  double precision :: u(1:nx+1), un1(1:nx+1), un2(1:nx+1)
  double precision :: tauair(1:nx+1)    ! tauair
  double precision :: b(1:nx+1)         ! b vector
  double precision :: zeta(0:nx+1), eta(0:nx+1), sigma(0:nx+1), Cw(1:nx+1), Cb(1:nx+1)
  double precision :: F_uk1(1:nx+1), R_uk1(1:nx+1) ! could use F for R
  double precision :: meanvalue, time1, time2, timecrap
  double precision :: L2norm, gamma_nl, nl_target, nbhr

  out_step = 0
  sigma    = 0d0 ! initial stresses are zero
  Nfail    = 0
  meanN    = 0

!------------------------------------------------------------------------
!     Input by user
!------------------------------------------------------------------------

  linear_drag    = .false.
  linear_viscous = .false. ! linear viscous instead of viscous-plastic
  constant_wind  = .true. ! T: 10m/s, F: spat and temp varying winds
  rep_closure    = .true. ! replacement closure (see Kreysher et al. 2000)
  restart        = .false.
  regularization = 'tanh' ! tanh, Kreyscher, capping (Hibler)
  adv_scheme     = 'upwind' ! upwind, upwindRK2
  oceanSIM       = .false. ! for shallow water model
  implicitDrag   = .true. ! for uwater mom eq.
  Asselin        = .true. ! Asselin filter for uw and etaw
  DiagStress     = .false. ! diagnostic for ice-ocean stress at i=idiag
  idiag          = 100
  Agamma         = 1d-02 ! Asselin filter parameter

  solver     = 2        ! 1: Picard+SOR, 2: JFNK, 3: EVP, 4: EVP*
  IMEX       = 0       ! 0: no IMEX, 1: Jdu=-F(IMEX), 2: J(IMEX)du=-F(IMEX) 
  BDF2       = 0       ! 0: standard, 1: Backward difference formula (2nd order)
  
  Deltat     = 600d0   ! time step [s]
  nstep      = 288     ! lenght of the run in nb of time steps
  Nmax_OL    = 200

  T = 0.36d0*Deltat ! elast. damping time scale (Deltate < T < Deltat)
  N_sub = 900
  Deltate    = Deltat / (N_sub*1d0) ! for EVP solver
  
  omega      = 1.5d0    ! relax parameter for SOR
  tol_SOR    = 1d-10    ! tol for SOR solver
  maxiteSOR  = 10000     ! max nb of ite for SOR
  iteSOR_pre = 10       ! nb of iterations for the SOR precond
  maxiteGMRES= 50      ! max nb of ite for GMRES
  gamma_nl = 1d-03
  dropini  = 1.5d0        ! defines initial drop in L2norm before gamma = 0.01
  small1   = 1d-10      ! to have a continuously diff water drag term
  small2   = 1d-22      ! to have a continuously diff rheology term
  smallA   = 1d-03      ! for num stab of Atw and Ata (in zones with ~no ice)

  expnb      = 1
  expres     = 2
  ts_res     = 50 ! time level of restart (!!! watchout for Deltat !!!)
  out_step(1)= 1   

!------------------------------------------------------------------------ 
! verify choice of solver and options
!------------------------------------------------------------------------ 

  if (solver .eq. 3 .or. solver .eq. 4) then
    if (IMEX .ne. 0 .and. BDF2 .ne. 0) then
      print *, 'set IMEX=0 and BDF2=0'
      stop
    endif
  endif
  if (BDF2 .eq. 1 .and. adv_scheme .ne. 'upwindRK2') then
      print *, 'set adv_scheme = upwindRK2'
      stop
  endif
  
!------------------------------------------------------------------------ 
!     Set first time level depending on restart specifications                
!------------------------------------------------------------------------

  if (restart) then
     tsini = ts_res + 1
  else
     tsini = 1
  endif
  
  tsfin = tsini - 1 + nstep
  
!------------------------------------------------------------------------
!     Define a flag for the precond (T) or solver (F)
!------------------------------------------------------------------------

  p_flag = .true.
  if (solver .eq. 1) p_flag = .false.

!------------------------------------------------------------------------
!     Define Deltax and check CFL based on input by user
!------------------------------------------------------------------------

  if ( nx .eq. 100 ) then 
     Deltax   =  20d03  ! grid size [m], the domain is always 2000 km 
  elseif  ( nx .eq. 200 ) then
     Deltax   =  10d03            
  elseif  ( nx .eq. 400 ) then
     Deltax   =  5d03            
  else
     print *,  'Wrong grid size dimenion', nx
     STOP
  endif

  Deltax2 = Deltax ** 2
  DtoverDx = Deltat / Deltax

  if ( 1d0*Deltat .gt. Deltax ) then
     print *, 'CFL condition is not respected'
     stop
  endif

!------------------------------------------------------------------------
!     Define constants
!------------------------------------------------------------------------

  C          = 20d0         ! ice strength parameter (watchout no A for now)
  Pstar      = 27.5d03      ! ice compression strength parameter
  e          = 2d0          ! ratio long to short axis of ellipse
  e_2        = 1/(e**2d0)   !
  alpha      = sqrt(1d0 + e_2)
  alpha2     = 1d0 + e_2
  kt         = 0d0          ! T = kt * P (1.0 in Konig and Holland, 2010)

  Cdair      = 1.2d-03      ! air-ice drag coeffient 
  Cdairw     = 1.2d-03      ! air-water drag coeffient 
  Cdwater    = 5.5d-03      ! water-ice drag coeffient
  rhoair     = 1.3d0        ! air density
  rho        = 900d0        ! ice density
  rhowater   = 1026d0       ! water density
  ge         = 9.8d0        ! Earth's gravitional acceleration
  Hw         = 2.5d0          ! mean water depth (for shallow water model)
  bw         = 0.0005d0/Hw  ! friction term for the uw momentum eq (always implicit).

  Cda        = rhoair   * Cdair
  Cdw        = rhowater * Cdwater

!------------------------------------------------------------------------
!     initial conditions
!------------------------------------------------------------------------

  call ini_get (u, restart, expres, ts_res)
  un1=u
  if ( BDF2 .eq. 1 ) un2 = un1 ! BDF2 needs u at 3 time levels
  tauair = 0d0 ! initialization (watchout for restart)
  nbhr = 0d0
  fgmres_per_ts = 0
  
  do ts = tsini, tsfin
     
     nbhr = nbhr + Deltat / 3600d0
     print *, 'time level, cumulative time (h) =', ts, nbhr
     
     
     call cpu_time(timecrap)
     call cpu_time(time1)

     if ( BDF2 .eq. 1 ) un2 = un1 ! BDF2 needs u at 3 time levels
                                  ! Attention not initialized the 1st time level.

!------------------------------------------------------------------------
!     update previous time level solutions
!------------------------------------------------------------------------

     un1=u
     hn1=h
     An1=A
     if (oceanSIM) then 
       uwn2   = uwn1 
       uwn1   = uw
       etawn2 = etawn1
       etawn1 = etaw
     endif
   
     if (IMEX .eq. 0) call ice_strength (hn1, An1) ! standard approach no IMEX 
     
!------- get wind forcing (independent of u) -----------------------------

     call wind_forcing (tauair, ts)
     
!------- Solves NL mom eqn at specific time step with solver1, 2 or 3
!        F(u) = A(u)u - b(u) = 0, u is the solution vector
!------- Begining of outer loop (OL) or Newton iterations ----------------
  
     if (solver .eq. 1 .or. solver .eq. 2 ) then ! implicit

	if (solver .eq. 2 ) call calc_scaling (An1) ! scaling=1 for other solvers
     
     do k = 1, Nmax_OL 
        
        if (IMEX .gt. 0) then ! IMEX method 1 or 2
	  call advection (un1, u, hn1, An1, h, A) ! advection scheme for tracers
	  call ice_strength (h, A) ! Pp_half is Pp/2 where Pp is the ice strength (Tp_half: tensile strength)
	endif
	
	call viscouscoefficient (u, zeta, eta) ! u is u^k-1
	call Cw_coefficient (u, Cw, Cb)            ! u is u^k-1
	call calc_R (u, zeta, eta, Cw, Cb, tauair, R_uk1)
	call Fu (u, un1, un2, h, R_uk1, F_uk1) 

        L2norm = sqrt(DOT_PRODUCT(F_uk1,F_uk1))
        
        if (k .eq. 1) then
	  nl_target = gamma_nl*L2norm
!	  call output_ini_L2norm(ts,L2norm,expnb)
	endif

	if (L2norm .lt. nl_target .or. L2norm .lt. 1d-08) exit

        if (solver .eq. 1) then
           print *, 'L2-norm after k ite=', ts, k-1, L2norm
           call bvect(tauair, un1, Cw, b)
           call SOR (b, u, h, A, zeta, eta, Cw, Cb, p_flag, ts)
!           call SOR_A (b, u, zeta, eta, Cw, k, ts)
        elseif (solver .eq. 2) then
           call prepFGMRES_NK(u, h, A, F_uk1, zeta, eta, Cw, Cb, un1, un2, tauair, &
                              L2norm, k, ts, fgmres_its)
!           call SOR_J(u, F_uk1, zeta, eta, Cw, upts, tauair, k, ts)
        endif
	fgmres_per_ts = fgmres_per_ts + fgmres_its
        if (k .eq. Nmax_OL) Nfail = Nfail + 1

     enddo
     meanN = meanN + k-1
!     call output_nb_ite (ts, k ,fgmres_per_ts, expnb)

     elseif (solver .eq. 3 .or. solver .eq. 4) then ! explicit (EVP)
     
      call viscouscoefficient (u, zeta, eta) ! u is u^k-1
      call Cw_coefficient (u, Cw, Cb)            ! u is u^k-1
      call EVP2solver(tauair, u, ts, solver)
     
     endif

!------------------------------------------------------------------------
!     Diagnostic of ice-ocean stress
!------------------------------------------------------------------------     
     if (DiagStress) then
      call Cw_coefficient (u, Cw, Cb)
      call calc_diag_stress (idiag, u, Cw, tauair)
     endif
!------------------------------------------------------------------------       
     
      call cpu_time(time2)
      print *, 'cpu time = ', time2-time1

     if (IMEX .eq. 0) call advection (un1, u, hn1, An1, h, A)  ! standard approach no IMEX

!     call meantracer(h,meanvalue)

!------------------------------------------------------------------------
!     Shallow water model
!------------------------------------------------------------------------
  
    if (oceanSIM) then
       
       call advect_etaw (etaw)
       call Cw_coefficient (u, Cw, Cb)
       if (IMEX .eq. 0) then
        call momentum_uw (idiag, tauair, Cdair, Cw, An1, u)
       elseif (IMEX .gt. 0) then
        call momentum_uw (idiag, tauair, Cdair, Cw, A, u)
       endif
       if (Asselin) then
	call Asselin_filter (etaw, uw)
       endif
    endif

!------------------------------------------------------------------------
!     output results
!------------------------------------------------------------------------

     if (ts .eq. out_step(1) .or. ts .eq. out_step(2) .or. &
         ts .eq. out_step(3) .or. ts .eq. out_step(4) .or. &
         ts .eq. out_step(5)) then
        print *, 'outputting results'
        call output_results(ts, expnb, solver, u, zeta, eta)
        call output_file(e, gamma_nl, solver, expnb)
     endif

!------------------------------------------------------------------------
!     calculate diagnostics            
!------------------------------------------------------------------------

!     call check_neg_vel(u)
     call minmaxtracer(h,1,ts)
     call minmaxtracer(A,2,ts)
     call minmaxtracer(u,3,ts)
     if (oceanSIM) then
      call minmaxtracer(etaw,4,ts)
      call minmaxtracer(uw,5,ts)
      if (DiagStress) call output_diag_stress (ts, expnb, idiag)
     endif

  enddo
  
  if (solver .eq. 1) then
     print *, 'Nb failures, mean ite of Picard: ', Nfail, (meanN*1d0)/(nstep*1d0)
  elseif (solver .eq. 2) then
     print *, 'Nb failures, mean ite of JFNK: ', Nfail, (meanN*1d0)/(nstep*1d0)
     print *, 'mean nb of fgmres it per time level: ', fgmres_per_ts/(nstep*1d0)
  endif

  deallocate(etaw, etawn1, etawn2, uw, uwn1, uwn2) 
  if (oceanSIM) then
   deallocate(duwdt, gedetawdx, tauiw, tauaw, buw) 
  endif
  
end program ice
      

