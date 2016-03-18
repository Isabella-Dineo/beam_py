      subroutine addnoise(iseed,array,n,amp)

      implicit none

      integer n, i, iseed
      real array(n),amp
      real ran1, gasdev, val

c      val = ran1(iseed)
c     assuming a normalized array
      write(*,*)'Seed for noise:', iseed
      if (amp.gt.0) then
         do i=1, n
            val = gasdev(iseed)*amp
            array(i) = array(i) + val
         enddo
      end if
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
CU    USES ran1
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software -0'>AX3'2500%'i6X@)1.
