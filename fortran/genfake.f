	program genfake

	implicit none
        include 'useful.inc'
	
	real rhomax, period_arr(7),ran1,be_max
	integer i,j,ngen(7)

        integer iseed,ncomp,npatch
        real*8 al,be,hmin,hmax,hrange,pw_in,cw
        real freq,period,noise
        logical verb,same_again
        character pgdev*50
	
	data period_arr/0.09,.3,.5,.7,.9,1.25,2./
	data ngen/16,71,46,24,24,27,18/
c	data ngen/100,100,100,100,100,100,100/

c     Get the input parameters
        call getparams(al,be,hmin,hmax,hrange,ncomp,npatch,
     +     cw,noise,period,freq,pgdev,iseed,verb)

        call pgbegin(0,pgdev,1,1)
        same_again=.false.

	do i=1,7
          if(period_arr(i).lt.0.1)then
	    hmin=950.
	    hmax=1000.
          else
	    hmin=20.
	    hmax=1000.
          endif
	  hrange=hmax-hmin
	  rhomax=3.*sqrt(pi*(hmax)/2./period_arr(i)/light)*180./pi
	  be_max=rhomax+3.*sqrt(rhomax)
	  be_max=rhomax
	  do j=1,ngen(i)
c choose alpha from 10 to 170
 5	    al=ran1(iseed)*180.
            if(al.lt.10. .or. al.gt.170.)goto 5
c choose beta from -be_max to +be_max
 10	    be=(2*ran1(iseed)-1.0)*be_max
	    if(abs(be).gt.be_max)goto 10

	    call pgadvance
            call beam_sub(al,be,hmin,hmax,hrange,ncomp,npatch,
     +          cw,noise,period_arr(i),freq,iseed,verb,same_again)

	  enddo
	enddo

	call pgend

	end
