      subroutine mapphi(al,be,phi,xp,yp)

      implicit none

      real*8 al, be, phi
      real*8 cosD, sinD, cosR, R, cosgamma, gamma, xp, yp
      real*8 acosD, asinD

      cosR = cosD(al+be)*cosD(al)
     +     +sinD(al+be)*sinD(al)*cosD(phi)
      call correct(cosR)
      R=acosD(cosR)
c     problems with precision for 180 deg
      if (int(r*100.0).eq.18000) R=int(R*100.0)/100.0
      if ((R.ne.0.0).and.(R.ne.180.0).and.(al.gt.0.0)) then
         cosgamma=(cosD(al+be)-cosD(al)*cosR)
     +        /(sinD(al)*sinD(R))
      else
         cosgamma=0.0
      endif
      call correct(cosgamma)
      gamma=acosD(cosgamma)
      xp=R*sinD(gamma)
      if (phi.gt.0) xp=-xp
      yp=-R*cosgamma 

      end

