
      subroutine prepFGMRES_NK(uk1, F_uk1, zeta, eta, Cw, upts, tauair, &
                               L2norm, k, ts, precond, fgmres_its)
        use size
        use numerical
        
      implicit none

      integer :: icode, iter, iout, i
      integer, intent(in) ::  k, ts, precond
      integer, intent(out) :: fgmres_its

      double precision, intent(inout) :: uk1(1:nx+1)
      double precision, intent(in)  :: F_uk1(1:nx+1), upts(1:nx+1)
      double precision, intent(in)  :: L2norm
      double precision, intent(in)  :: zeta(0:nx+1), eta(0:nx+1)
      double precision, intent(in)  :: Cw(1:nx+1)
      double precision, intent(in) :: tauair(1:nx+1) 
      double precision :: du(1:nx+1), rhs(1:nx+1)
      double precision :: vv(1:nx+1,img1), wk(1:nx+1,img)!, Funeg(1:nx+1)
      double precision :: wk1(1:nx+1), wk2(1:nx+1)
      double precision :: eps, gamma, epsilon, s

!------------------------------------------------------------------------
!     This routine solves J(u)du = -F(u) where u = u^k, du = du^k using the
!     Jacobian free Newton Krylov method. The Krylov method is the precon-
!     ditioned FGMRES method. The usefull references are:
!
!     Knoll and Keyes, J.of.Comput.Physics, 2004.
!     Eisenstat and Walker, SIAM J.Optimization, 1994.
!     Eisenstat and Walker, SIAM J.Sci.Comput., 1996.
!
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!     Making of the RHS vector : -F(uk1) and calculation of res norm
!------------------------------------------------------------------------

      rhs = -1d0*F_uk1 ! mult by -1 because we solve Jdu = -F(u)

!------------------------------------------------------------------------
!     Initial guess vector: because we solve for du and du is just a 
!     correction, we set the initial guess to zero
!------------------------------------------------------------------------

      du = 0d0

!------------------------------------------------------------------------
!     Choosing the forcing term (eta_e)
!------------------------------------------------------------------------
      
      call forcing_term (k,ts,L2norm,gamma)

!------------------------------------------------------------------------
!      Begining of FGMRES method    
!------------------------------------------------------------------------

      eps = gamma * L2norm ! setting the tolerance for fgmres

      iout   = 0    ! set  higher than 0 to have res(ite)

      icode = 0

 10   CONTINUE
      
      call fgmres (nx+1,img,rhs,du,iter,vv,wk,wk1,wk2, &
                   eps,maxiteGMRES,iout,icode,fgmres_its)

      IF ( icode == 1 ) THEN
!         CALL identity (wk1,wk2)
         if (precond .eq. 1) then
            CALL SOR (wk1, wk2, zeta, eta, Cw, .true., ts)
         elseif (precond .eq. 2) then
!            CALL EVP1Bprecond(wk1, wk2, zeta, eta, Cw, ts)
            CALL EVP2Dprecond(wk1, wk2, zeta, eta, Cw, ts)
         endif

         GOTO 10
      ELSEIF ( icode >= 2 ) THEN
         epsilon = 1d-07 ! approximates Jv below
         call JacfreeVec (wk1, wk2, F_uk1, uk1, upts, tauair, epsilon) 
         GOTO 10
      ENDIF

!------------------------------------------------------------------------
!      End of FGMRES method    
!------------------------------------------------------------------------

      if (fgmres_its .eq. maxiteGMRES) then
         print *,'WARNING: FGMRES has not converged'
         print*, 'Please check the precond relaxation param (wlsor or wsor).'
         stop
      endif

! icode = 0 means that fgmres has finished and sol contains the app. solution

!------------------------------------------------------------------------
!      line search method: add a*du to u (where a = 0.25, 0.5 or 1.0)
!------------------------------------------------------------------------

!         call linesearch(sol, x, res)
!         print *, 'res after linesearch = ', res, eta_e
!	 call output_u_and_du ( ts, k, uk1, du )
	 call calc_s( uk1, du, s )
         uk1 = uk1 + du ! u^k+1 = u^k + du^k

         return
       end subroutine prepFGMRES_NK
      
   subroutine forcing_term(k,ts,L2norm,gamma)
  
      use numerical
      implicit none

      integer, intent(in) :: k, ts

      double precision, intent(in) :: L2norm
      double precision, save :: L2normk_1, L2norm_t
      double precision :: gamma_ini, phi_e, alp_e
      double precision, intent(out) :: gamma

      gamma_ini = 0.99d0
      phi_e     = 1d0
      alp_e     = 1d0 !      alp_e = (1d0 + 5d0**0.5d0)/2d0 !2d0 

      if (k .eq. 1) then

         gamma = gamma_ini
         L2norm_t = L2norm / dropini ! t stands for transition

      elseif (k .gt. 200) then

         gamma = gamma_ini

      else

         if (L2norm .gt. L2norm_t) then

            gamma = gamma_ini
      
         else

            gamma = phi_e * (L2norm/L2normk_1)**alp_e ! Eisenstat, 1996,eq2.6 
            gamma = min(gamma_ini,gamma)
            gamma = max(0.01d0,gamma)

         endif

      endif

      L2normk_1 = L2norm

      if (ts .le. 10) gamma = gamma_ini

    end subroutine forcing_term


   subroutine calc_s(uk1, du, s)
      use size
      use numerical
      implicit none

      integer :: i, imaxdu

      double precision, intent(in) :: uk1(1:nx+1), du(1:nx+1)
      double precision, intent(out) :: s
      double precision :: aa, temp, maxdu

      aa = 0.25d0
      temp=100000d0
      maxdu = 0d0
      
      do i = 2, nx
	if ( abs(aa*uk1(i)/du(i)) .lt. temp) temp =abs(aa*uk1(i)/du(i))
	if ( abs(du(i)) .gt. maxdu ) then
	  maxdu = abs(du(i))
	  imaxdu = i
	endif
      enddo
  
      s = min(1d0,temp)

    end subroutine calc_s

