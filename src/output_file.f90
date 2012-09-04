subroutine output_file(e, gamma_nl, solver, precond, expnb)
 
  use size
  use rheology
  use resolution
  use numerical
  use EVP_const
  use option

  implicit none

  character filename*30

  integer, intent(in) :: solver, precond, expnb
  double precision, intent(in) :: e, gamma_nl
  
  write (filename, '("output/info.",i2.2)') expnb
  open (10, file = filename, status = 'unknown')
  
  write(10,*) ('PHYSICAL PARAMETERS')
  write(10,*) ('')
  write(10,*) ('C =')
  write(10,10) (C)
  write(10,*) ('Pstar =')
  write(10,10) (Pstar)
  write(10,*) ('e =')
  write(10,10) (e)
  write(10,*) ('zetamin =')
  write(10,10) (zetamin)
  write(10,*) ('zetamax parameter =')
  write(10,10) (zmax_par)

  write(10,*) ('')
  write(10,*) ('DEFINITION OF THE RUN')
  write(10,*) ('')

  write(10,*) ('dx (m) =')
  write(10,10) (Deltax)
  write(10,*) ('dt (s) =')
  write(10,10) (Deltat)
  write(10,*) ('linear water drag =')
  write(10,*) (linear_drag)
  write(10,*) ('linear viscous =')
  write(10,*) (linear_viscous)
  write(10,*) ('constant wind =')
  write(10,*) (constant_wind)

  write(10,*) ('')
  write(10,*) ('SOLVER AND NUMERICS')
  write(10,*) ('')

  if (implicit_solv) then
     write(10,*) ('SOLVER = ')
     write(10,11) (solver)
     write(10,*) ('precond = ')
     write(10,11) (precond)
     if (precond .eq. 1) then
        write(10,*) ('omega (sor) =')
        write(10,10) (omega)
        write(10,*) ('nb of sor ite =')
        write(10,12) (iteSOR_pre)
     elseif (precond .eq. 2) then
        write(10,*) ('EVP precond')
        write(10,*) ('dte (s) =')
        write(10,10) (Deltate)
        write(10,*) ('nb of sub cycles =')
        write(10,12) (N_sub)
     endif
     write(10,*) ('max nb of gmres ite = ')
     write(10,12) (maxiteGMRES)
     write(10,*) ('gamma_nl (tol) =')
     write(10,10) (gamma_nl)
     
  elseif (.not. implicit_solv) then
     write(10,*) ('SOLVER = EVP')
     write(10,*) ('dte (s) =')
     write(10,10) (Deltate)
     write(10,*) ('T (s) =')
     write(10,10) (T)
  endif

  close(10)

10 format (f24.10)
11 format (i1.1)
12 format (i3.3)

  return
end subroutine output_file


