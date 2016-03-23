	program genpop

	implicit none
        include 'useful.inc'

	integer nmax
	parameter (nmax=5000)
	
	real rhomax, ran1
	real period(nmax),age(nmax),edot(nmax)
	real area1,area2,area,iptot
	real ymin,ymax,hist(nmax),xhist(nmax)
	integer i,j,ntot,ngen,tgood
	real single,double,multiple,nodet

        integer iseed,ncomp,npatch
        real*8 al,be,zeta,hmin,hmax,cw,decay
	real*8 al_ip,be_ip,alpha(nmax),beta(nmax)
	real dm
	real pms(nmax),pdot15(nmax)
	real al2(nmax),be2(nmax)
        real freq,noise,fwhm_save,wmax_save
        logical same_again,flatpol,dobeam,doprof,doab
        character pgdev*50,iplab*2,key_save*10,outfile*20

c     Get the input parameters
        call getparams(al,be,hmin,hmax,ncomp,npatch,cw,noise, period(1)
     +    ,freq,pgdev,iseed,verb,flatpol,dobeam,doprof,doab,dm,
     +    outfile)

c	open(unit=10,file='highedot.list',form='formatted',status='old')
	open(unit=10,file='dunc.list',form='formatted',status='old')
	i=1
c 5      read(10,'(15x,f8.0,f18.0,f18.0)',end=10)period(i),age(i),edot(i)
 5      read(10,*,end=10) pms(i),pdot15(i),alpha(i),beta(i)
	i=i+1
	goto 5

 10	ntot=i-1
	ngen=1
	if(ngen*ntot .gt. nmax)then
	  write(*,'('' Too many trials !! '')')
	  stop
	endif
	iptot=0.0
	single=0.0
	double=0.0
	multiple=0.0
	nodet=0.0
	
        call pgbegin(0,pgdev,1,1)
        same_again=.false.

	do i=1,ntot

	  write(*,'('' Loop '',i5,'' of '',i5)')i,ntot*ngen

c set hmin and hmax according to the Edot
c and also the number of components and patches
	  edot(i)=31+log10(3.95*pdot15(i)*(pms(i)/1000.)**(-3))
	  age(i)=15.8e+06*pms(i)/1000/pdot15(i)
	  write(*,*) 'Paras:',edot(i),age(i),alpha(i),beta(i),pms(i)
	  if (edot(i) .lt. 35.0)then
	    hmin=50.
	    hmax=500.
	    ncomp=4
	    npatch=4
	  else
	    hmin=450.
	    hmax=500.
	    ncomp=2
	    npatch=5
	  endif
	  period(i)=pms(i)/1000.
	  rhomax=3.*sqrt(pi*(hmax)/2./period(i)/light)*180./pi
	  decay=exp((10**age(i))/-70.e6)
	  al=alpha(i)
	  be=beta(i)
	  call pgadvance
	  call beam_sub(al,be,hmin,hmax,ncomp,npatch,cw,noise,
     +          period(i),freq,iseed,same_again,flatpol,
     +          dobeam,doprof,doab,dm)

	  if(key .eq. 'S')single=single+1.
	  if(key .eq. 'D')double=double+1.
	  if(key .eq. 'M')multiple=multiple+1.
	enddo
	call pgend
	
	end
