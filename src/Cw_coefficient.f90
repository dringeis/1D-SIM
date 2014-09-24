subroutine Cw_coefficient (utp, Cw, Cb)
  use size
  use forcing
  use option  
  use global_var
  
  implicit none

  integer :: i

  double precision, intent(in)  :: utp(1:nx+1)
  double precision, intent(out) :: Cw(1:nx+1), Cb(1:nx+1) ! defined at u location
  double precision :: utypical, A_at_u, h_at_u, bathy_at_u, umin
  double precision :: Cbfactor, k1, k2, hc, CC ! similar to C in ice_strength calc...in fact we set CC=C=20d0
  
  utypical = 0.1d0
  CC = 20d0
  k1 = 10d0
  k2 = 10d0
  umin=1d-03
  
  Cw(1)    = 0d0
  Cw(nx+1) = 0d0
  Cb(1)    = 0d0
  Cb(nx+1) = 0d0

  if (linear_drag) then
     
     Cw(2:nx) = Cdw*utypical     ! linear water drag
     
  else
     
     do i = 2, nx
        
        Cw(i) = Cdw*sqrt(utp(i)**2d0 + small1)
        
!        bathy_at_u = ( bathy(i-1) + bathy(i) ) / 2d0
	bathy_at_u = min(bathy(i-1), bathy(i))
	
        if (bathy_at_u .gt. 30d0) then ! too deep for bottom drag
	  Cb(i) = 0d0
	else
	  hc=bathy_at_u/k1
	  h_at_u = max( h(i-1), h(i) )
	  
	  if ( h_at_u .gt. hc ) then
	    A_at_u = max( A(i-1), A(i) )
	    Cbfactor=k2/(abs(utp(i))+umin)
	    Cb(i) = Cbfactor * (h_at_u - hc) * dexp(-CC * ( 1d0 - A_at_u ))
	  else
	    Cb(i) = 0d0
	  endif
	 endif

     enddo
     
  endif

  return
end subroutine Cw_coefficient



