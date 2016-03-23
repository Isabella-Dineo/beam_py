      integer nz, np, compmax
      real pi,light
      parameter (nz=2000, np=1024, compmax=100)
      parameter (pi=3.1415926545, light=3.e5)

      real*8 pcd(compmax),hcomp(compmax),pw(compmax)
      real*8 new_hcomp(compmax),rho(compmax)

      real profi(np),profq(np),profu(np),xprof(np),profl(np),profv(np)
      real profpa(np),epa(np),xlos(np),ylos(np),thetalos(np)
      real fwhm,widthmax
      real profxmin,profxmax,profymin,profymax
      character key*10
      logical verb

      common /c1/ pcd,hcomp,pw,new_hcomp,rho
      common /c2/ profi,profq,profu,xprof,profl,profv,
     +            profpa,epa,xlos,ylos,thetalos,
     +            profxmin,profxmax,profymin,profymax,
     +            fwhm,widthmax,verb,key
