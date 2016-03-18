	subroutine rvm(alpha,beta,phi0,psi0,xprof,pa,np)

	implicit none
	real pi
	parameter (pi=3.1415926545)
	integer np,i
	real*8 alpha,beta,phi0,psi0,zeta,phi
	real xprof(np),pa(np),num,denom

c convert to radians
	alpha=alpha*pi/180.
	beta=beta*pi/180.

c psi convention problem
c	beta=-beta
c	alpha=pi-alpha
c

	zeta=beta+alpha
	do i=1,np
	  phi=xprof(i)*pi/180.
	  num=sin(alpha)*sin(phi-phi0)
	  denom=sin(zeta)*cos(alpha)-cos(zeta)*sin(alpha)
     +          *cos(phi-phi0)
	  pa(i)=psi0 + atan(num/denom)
	  if(pa(i).lt.-pi/2.)pa(i)=pa(i)+pi
	  if(pa(i).gt.pi/2)pa(i)=pa(i)-pi
c convert back to degrees
          pa(i)=pa(i)*180./pi
	enddo

	return
	end
