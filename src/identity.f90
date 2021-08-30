
subroutine identity ( rhs, wk2 ) ! jfl first one is the
                                 ! solution after preconditioning

  use size

!--------------------------------------------------------------------------
! Preconditioner: calculates M x = rhs where rhs is wk1. x is the initial
! guess. x is set to 0 here. The solution x is then put in wk2.
!--------------------------------------------------------------------------

  implicit none

  double precision, intent(in)  :: rhs(1:nx+1)
  double precision, intent(out) :: wk2(1:nx+1)

  wk2 = rhs

  return
end subroutine identity






