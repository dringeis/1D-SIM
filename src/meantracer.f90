subroutine meantracer(var, meanvalue)
  use size
  
  implicit none
  
  integer :: i

  double precision, intent(in) :: var(0:nx+1)
  double precision, intent(out) :: meanvalue
  external fluxx
  
  meanvalue = 0d0
  
  do i = 1, nx

     meanvalue = meanvalue + var(i)

  enddo

  meanvalue = meanvalue / (nx*1d0)
  
  return
end subroutine meantracer









