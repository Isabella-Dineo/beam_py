      program beam

      implicit none

      include 'useful.inc'

      integer iseed,ncomp,npatch,i,imjd,mjds
      real*8 al,be,hmin,hmax,cw,mjdf
      real freq,period,noise,dm,flux
      logical same_again,flatpol,dobeam,doprof,doab
      character pgdev*50,outfile*20, template*14,outfits*24

c     Get the input parameters
      call getparams(al,be,hmin,hmax,ncomp,npatch,cw,noise,
     +     period,freq,pgdev,iseed,verb,flatpol,dobeam,doprof,doab,dm,
     +     outfile,flux)

      call pgbegin(0,pgdev,1,1)
      same_again=.false.
      call beam_sub(al,be,hmin,hmax,ncomp,npatch,cw,noise, period,freq
     $     ,iseed,same_again,flatpol,dobeam,doprof,doab,dm,flux)
      call pgend
      open(unit=77,file=outfile, form='formatted', status='unknown')
      do i=1,np
         write(77,*)xprof(i),i*period*1000/np,profi(i),profq(i),profu(i)
     $        ,profv(i)
      end do
      write(*,*)
      close(77)
c     Experimental - Write out fits file
      template='testtemp.txt'
      outfits='psrsim.fits'
      imjd=54867 ! 5 Feb 2009
      mjds=0
      mjdf=0.0
      call beamfits(period,freq,template,outfits,imjd,mjds,mjdf)
      end
