subroutine Cw_coefficient (utp, Cw, Cb)
  use size
  use forcing
  use option  
  use global_var
  
  implicit none

  integer :: i

  double precision, intent(in)  :: utp(1:nx+1)
  double precision, intent(out) :: Cw(1:nx+1), Cb(1:nx+1)
  double precision :: utypical, A_at_u, h_at_u, umin
  double precision :: Cbfactor, KK, CC ! similar to C in ice_strength calc...in fact we set CC=C=20d0
  
  utypical = 0.1d0
  CC = 20d0
  KK = 10d0
  umin=1d-02
  
  Cw(1)    = 0d0
  Cw(nx+1) = 0d0
  Cb(1)    = 0d0
  Cb(nx+1) = 0d0

  if (linear_drag) then
     
     Cw(2:nx) = Cdw*utypical     ! linear water drag
     
  else
     
     do i = 2, nx
        
        Cw(i) = Cdw*sqrt(utp(i)**2d0 + small1)
        if (bathy(i) .gt. 20d0) then ! too deep for bottom drag
	  Cb(i) = 0d0
	else
	  A_at_u = ( A(i-1) + A(i) ) / 2d0
	  h_at_u = ( h(i-1) + h(i) ) / 2d0
	  Cbfactor=KK/(abs(utp(i))+umin)
	  Cb(i) = Cbfactor * h_at_u * dexp(-CC * ( 1d0 - A_at_u ) )
	  !Cb(i) = Cbfactor
	endif

     enddo
     
  endif

  return
end subroutine Cw_coefficient



