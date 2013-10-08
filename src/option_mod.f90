MODULE option
!
  IMPLICIT NONE
  logical :: linear_drag, linear_viscous, constant_wind
  logical :: rep_closure
  integer :: IMEX, BDF2
  character(LEN=20) :: adv_scheme, regularization

END MODULE option
