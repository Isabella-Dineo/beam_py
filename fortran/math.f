      real*8 function cosD(arg)
      
      real*8 arg, degtorad
      parameter (degtorad=0.01745329)  ! pi/

      cosD=cos(arg*degtorad)

      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function sinD(arg)
      
      real*8 arg, degtorad
      parameter (degtorad=0.01745329)  ! pi/

      sinD=sin(arg*degtorad)

      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function asinD(arg)
      
      real*8 arg, degtorad
      parameter (degtorad=0.01745329)  ! pi/

      asinD=asin(arg)/degtorad

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function acosD(arg)
      
      real*8 arg, degtorad
      parameter (degtorad=0.01745329)  ! pi/

      acosD=acos(arg)/degtorad

      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function atanD(arg)
      
      real*8 arg, degtorad
      parameter (degtorad=0.01745329)  ! pi/

      atanD=atan(arg)/degtorad

      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine correct(x)

c     coorect for round-off errors

      implicit none

      real*8 x, sign, y, tol

      tol=1.0e-7

      sign=1.0
      if (x.lt.0.0) sign=-1.0

      y=abs(x)
      if (abs(y-1.0).lt.tol) y=1.0
      if (abs(y-0.0).lt.tol) y=0.0

      x=y*sign
      
      end

