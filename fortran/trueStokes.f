      program trueStokes
c     generate a profile from beam, then create a number of frequency
c     channels with the correct dispersion/Faraday rotation
      implicit none

      include 'useful.inc'

      integer iseed,ncomp,npatch,i
      real*8 al,be,hmin,hmax,cw
      real freq,period,noise
      real dm,rm
      logical same_again,flatpol,dobeam,doprof,doab
      character pgdev*50

c     Get the input parameters
      call getparams(al,be,hmin,hmax,ncomp,npatch,cw,noise, period,freq
     $     ,pgdev,iseed,verb,flatpol,dobeam,doprof,doab)

      call pgbegin(0,pgdev,1,1)
      same_again=.false.     

      call beam_sub(al,be,hmin,hmax,ncomp,npatch,cw,noise, period,freq
     $     ,iseed,same_again,flatpol,dobeam,doprof,doab)
      call pgend

c     Ask for dm and rm [later pick them up from catalogue]
      write(*,*) 'Give the dm,rm in usual units:'
      read (*,*) dm,rm
      write(*,*) dm,rm
      do i=1,np
         write(*,*) i,profi(i),profq(i),profu(i),profv(i)
      end do

      end
      
