      subroutine plotprof(al,be,xopt,yopt)

      implicit none

      include 'useful.inc'

      real*8 al,be
      real ymin,ymax
      character angle*10,ach*20,bch*20,kch*20
      character xopt*10,yopt*10

      write(angle,'(f4.0,'' '',f4.0)') al,be
      write(ach,'('' Alpha = '',f4.0)')al
      write(bch,'('' Beta  = '',f4.0)')be
      write(kch,'('' Class = '',a1)')key

c plot the profile
      ymax=1.1*profymax
      ymin=0.0
      if(profymin.lt.0.)ymin=1.1*profymin
c      profxmin=-180.
c      profxmax=180.
      call pgswin(profxmin,profxmax,ymin,ymax)
      call pgsci(1)
      call pgslw(1)
         call pgbox(xopt,0.0,0,yopt,0.0,0)      
         call pgslw(1)
         call pglabel('Phase (deg)','Intensity',' ')
         call pgslw(3)
         call pgsci(1)
         call pgline(np,xprof,profi)
         call pgsci(2)
         call pgline(np,xprof,profl)
         call pgsci(4)
         call pgline(np,xprof,profv)
	 call pgsvp(0.08,0.49,0.1,0.37)
	 call pgswin(0.,1.,0.,1.)
         call pgsci(1)
         call pgtext(0.1,0.4,ach)
         call pgtext(0.1,0.3,bch)
         call pgtext(0.1,0.2,kch)
      
      return
      end

