	subroutine getlos(al,be,xlos,ylos,thetalos,np)

	implicit none
	integer i,np
	real*8 al,be,xp,yp,phi,step
	real xlos(np),ylos(np),thetalos(np)

	phi=-180.
	step=360./float(np)
	do i=1,np
	  phi=phi+step
	  call mapphi(al, be, phi, xp, yp)
	  xlos(i)=xp
	  ylos(i)=yp
          thetalos(i)=atan2(yp,xp)*180./3.14159265 - 90.0
	enddo

	return
	end
