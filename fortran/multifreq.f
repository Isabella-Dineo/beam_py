      program multifreq

      implicit none

      include 'useful.inc'

      integer iseed,ncomp,npatch,ff
      real*8 al,be,hmin,hmax,cw
      real freq,period,noise
      real freqs(10),bly,try
      logical same_again,flatpol,dobeam,doprof
      character pgdev*50,xopt*10,yopt*10

c     Get the input parameters
      call getparams(al,be,hmin,hmax,ncomp,npatch,cw,noise,
     +     period,freq,pgdev,iseed,verb,flatpol,dobeam,doprof)

      call pgbegin(0,pgdev,1,1)

      dobeam=.true.
      doprof=.false.
      xopt='bc'
      yopt='bc'
      do ff=1,10
        if(ff.eq.1)then
          freqs(1)=0.1
          same_again=.false.
        else
          freqs(ff)=freqs(ff-1)*1.7
          same_again=.true.
        endif
        call beam_sub(al,be,hmin,hmax,ncomp,npatch,cw,noise,
     +     period,freqs(ff),iseed,same_again,flatpol,dobeam,doprof)
        try=1.05-0.09*real(ff)
        bly=try-0.09
        if(ff.eq.10)xopt='bcnst'
        call pgsvp(0.57,0.99,bly,try)
        call plotprof(al,be,xopt,yopt)

	write(*,*)profxmin,profxmax,profymin,profymax
      enddo

      call pgend

      end
