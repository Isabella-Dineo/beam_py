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
	real*8 al_ip,be_ip
	real al2(nmax),be2(nmax)
        real freq,noise,fwhm_save,wmax_save
        logical same_again,flatpol,dobeam,doprof,doab
        character pgdev*50,iplab*2,key_save*10

c     Get the input parameters
        call getparams(al,be,hmin,hmax,ncomp,npatch,cw,noise,
     +     period(1),freq,pgdev,iseed,verb,flatpol,dobeam,doprof,doab)

c	open(unit=10,file='highedot.list',form='formatted',status='old')
	open(unit=10,file='alledot2.list',form='formatted',status='old')
	read(10,*)
	i=1
 5      read(10,'(15x,f8.0,f18.0,f18.0)',end=10)period(i),age(i),edot(i)
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
	  if (edot(i) .lt. 35.0)then
	    hmin=20.
	    hmax=1000.
	    ncomp=4
	    npatch=5
	  else
	    hmin=950.
	    hmax=1000.
	    ncomp=1
	    npatch=10
	  endif

	  rhomax=3.*sqrt(pi*(hmax)/2./period(i)/light)*180./pi
	  decay=exp((10**age(i))/-70.e6)

	  do j=1,ngen

c choose random birth alpha (cos distribution)
 20 	    al=acos(ran1(iseed))*180./pi
	    al=al*decay

c choose random birth zeta (cos distribution)
c and hence beta
c check for detection - if not regenerate
  	    zeta=acos(ran1(iseed))*180./pi
	    be=zeta-al
	    if(abs(be).gt.rhomax)goto 20

 	    call pgadvance
            call beam_sub(al,be,hmin,hmax,ncomp,npatch,cw,noise,
     +          period(i),freq,iseed,same_again,flatpol,
     +          dobeam,doprof,doab)

	    if(key .eq. 'S')single=single+1.
	    if(key .eq. 'D')double=double+1.
	    if(key .eq. 'M')multiple=multiple+1.
	    iplab(1:1)='Y'
	    fwhm_save=fwhm
	    key_save=key
	    wmax_save=widthmax
c main pulse not detected 
	    if(key .eq. '-')iplab(1:1)='N'

c check for interpulse detection
c note new alpha and beta (see p182 in Patrick's thesis)
c then make a beam for that pole
	    if(abs(180. - 2*al - be).lt.rhomax)then
	      al_ip = 180. - al
	      be_ip = be + 2.*al - 180.
 	      call pgadvance
              call beam_sub(al_ip,be_ip,hmin,hmax,ncomp,npatch,cw,noise,
     +            period(i),freq,iseed,same_again,flatpol,
     +            dobeam,doprof,doab)

	      if(key .eq. '-')then
                iplab(2:2)='N'
              else
                iplab(2:2)='Y'
	        fwhm_save=fwhm
	        key_save=key
	        wmax_save=widthmax
              endif
            else
              iplab(2:2)='N'
            endif

c if nothing has been detected try again
	    if(iplab.eq.'NN')then
	      nodet=nodet+1.0
              goto 20
	    endif
	    if(iplab.eq.'YY')then
	      iptot=iptot+1.0
	    endif

	    write(33,'(7f8.3,x,a,x,a)')period(i),age(i),edot(i),
     +        al,be,fwhm_save,wmax_save,iplab,key_save
	    al2(i)=al
	    be2(i)=be

	  enddo
	enddo

	write(*,'('' number generated '',i6)')ntot*ngen
	write(*,'('' non dets '',i6,'' = ''f6.2,'' %'')')nint(nodet),
     +     nodet/ntot/ngen*100.
	write(*,'('' singles '',i6,'' = ''f6.2,'' %'')')nint(single),
     +     single/ntot/ngen*100.
	write(*,'('' doubles '',i6,'' = ''f6.2,'' %'')')nint(double),
     +     double/ntot/ngen*100.
	write(*,'('' multiples '',i6,'' = ''f6.2,'' %'')')nint(multiple),
     +     multiple/ntot/ngen*100.
	write(*,'('' interpulses '',i6,'' = ''f6.2,'' %'')')nint(iptot),
     +     iptot/ntot/ngen*100.

	call dohist(al2,ntot*ngen,0.0,90.0,ymin,ymax,tgood,hist,xhist,20)

	call pgend

	end
