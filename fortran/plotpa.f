      subroutine plotpa(pa_rvm)

      implicit none
      include 'useful.inc'

      integer i,gb
      real ey1,ey2,pa_rvm(np)

c     -- PA profile
      call pgsvp(0.57,0.99,0.70,0.95)
      call pgswin(profxmin,profxmax,-90.0,90.0)
      call pgsci(1)
      call pgslw(3)
      call pgbox('BCST',0.0,0,'BCNST',0.0,0)      
      call pgslw(1)
      call pglabel('','PA (deg)','')

      gb=-1
      do i=1,np
	if(profpa(i).ne.0.0)then
          call pgpt(1,xprof(i),profpa(i),16)
	  ey1=profpa(i)-epa(i)
	  if(ey1.lt.-90)ey1=-90.0
	  ey2=profpa(i)+epa(i)
	  if(ey2.gt.+90)ey2=+90.0
	  call pgmove(xprof(i),ey1)
	  call pgdraw(xprof(i),ey2)
	  if(i.eq.gb+1)then
	    call pgmove(xprof(i),profpa(i))
	    if(abs(profpa(i)-profpa(gb)).le.90.)
     +        call pgdraw(xprof(gb),profpa(gb))
	  endif
	  gb=i
	endif
      enddo

      call pgsci(2)
      call pgline(np,xprof,pa_rvm)
c      call pgsci(3)
c      call pgline(np,xprof,thetalos)

      end
