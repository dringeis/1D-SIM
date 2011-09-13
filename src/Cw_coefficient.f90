subroutine Cw_coefficient (utp, Cw)
  use size
  use forcing
  use option
  
  implicit none

  integer :: i

  double precision, intent(in)  :: utp(1:nx+1)
  double precision, intent(out) :: Cw(1:nx+1)
  double precision :: utypical
  
  utypical = 0.1d0

  Cw(1)    = 0d0
  Cw(nx+1) = 0d0

  if (linear_drag) then
     
     Cw(2:nx) = Cdw*utypical     ! linear water drag
     
  else
     
     do i = 2, nx
        
        Cw(i) = Cdw*abs(utp(i))  ! nonlinear water drag
        
     enddo
     
  endif

  return
end subroutine Cw_coefficient



