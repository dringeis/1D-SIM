MODULE rheology
!
! C : ice strength parameter
! Pstar : ice compression strength parameter
! zetamin : minimum bulk viscosity 
! ell_2 : 1/ellipticity**2
! zmax_par : capping value of zeta 
! zmax_par_1 : 1/zmax_par
  IMPLICIT NONE
  DOUBLE PRECISION :: C, Pstar, alpha, alpha2
  DOUBLE PRECISION :: zetamin, e_2, zmax_par

END MODULE rheology

MODULE properties

  IMPLICIT NONE
  DOUBLE PRECISION :: rho

END MODULE properties

MODULE forcing

  IMPLICIT NONE
  DOUBLE PRECISION :: Cda, Cdw

END MODULE forcing

MODULE resolution

  IMPLICIT NONE
  DOUBLE PRECISION :: Deltax, Deltax2, Deltat, Deltate

END MODULE resolution

MODULE numerical

  IMPLICIT NONE
  INTEGER :: N_sub, maxiteSOR, maxiteGMRES, iteSOR_pre
  DOUBLE PRECISION :: T
  DOUBLE PRECISION :: omega, tol_SOR

END MODULE numerical

MODULE EVP_const

  IMPLICIT NONE
  DOUBLE PRECISION :: beta_1, mindeno, one_or_zero, Estar

END MODULE EVP_const
