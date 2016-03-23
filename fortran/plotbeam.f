      subroutine plotbeam(zp,nz,plotxmin,plotxmax)

      implicit none

      integer nz, nz4, i, j
      parameter (nz4=2000)
      real*8 zp(nz,nz)
      real*8 sinD, cosD, ang
      real zmin,zmax
      real matrix(6), fg, bg
      real x1,x2,y1,y2,plotxmin,plotxmax
      real xx1, yy1, zp4(nz,nz)

      if (nz4.lt.nz) stop 'plotbeam: array too small'

c transfer to real*4 and compute zmin and zmax
      zmin=0.0
      zmax=0.0
      do i=1,nz
         do j=1, nz
            zp4(i,j)=zp(i,j)
            if(zp4(i,j).gt.zmax)zmax=zp4(i,j)
         enddo
      enddo
 
c looks better on gray scale
      zmax=zmax/2.

      x1=-180.
      x2=180.
      y1=-180.
      y2=180.

c      x=m1+m2*i+m3*j
      matrix(1)=x1
      matrix(2)=(x2-x1)/float(nz)
      matrix(3)=0.0

c      y=m4+m5*i+m6*j
      matrix(4)=y1
      matrix(5)=0.0
      matrix(6)=(y2-y1)/float(nz)

      fg=zmax
      bg=zmin+(0.1*(zmax-zmin))

      call pgsci(1)
      call pgslw(1)
      call pgsvp(0.08,0.49,0.37,0.95)
      call pgswin(plotxmin,plotxmax,plotxmax,plotxmin)
      call pgbox('BCNST',0.0,0,'BCNST',0.0,0)      
      call pglabel('X (deg)','Y (deg)',' ')

      call pggray(zp4,nz4,nz4,1,nz,1,nz,fg,bg,matrix)

      call pgsls(2)
      do i=0,360,20
	ang=i
	xx1=500.*sinD(ang)
	yy1=500.*cosD(ang)
	call pgmove(0.0,0.0)
	call pgdraw(xx1,yy1)
      enddo
      call pgsls(1)

      end


