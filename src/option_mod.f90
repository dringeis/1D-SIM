MODULE option
!
  IMPLICIT NONE
  logical :: linear_drag, linear_viscous, constant_wind
  logical :: rep_closure, implicit_solv
  integer :: IMEX
  character(LEN=20) :: adv_scheme

END MODULE option
