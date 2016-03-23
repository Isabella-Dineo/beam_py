	program genfake

	implicit none
	real pi,light,hmax
	parameter (pi=3.1415926545, light=3.e5, hmax=1000.)
	
	real rhomax, period(7),ran1,alpha,beta,be_max
	integer i,j,ngen(7),iseed,i_alpha,i_beta
	character extras*50
	
c	data period/0.11,.3,.5,.7,.9,1.25,2./
	data period/0.09,.3,.5,.7,.9,1.25,2./
	data ngen/16,71,46,24,24,27,18/
c	data ngen/100,100,100,100,100,100,100/

	write(*,'('' Input seed >''$)')
	read(*,*)iseed
	if(iseed.gt.0)iseed=-iseed
	alpha=ran1(iseed)

	open(unit=10,file='fake.file',status='unknown')

	extras='-na0.1 -nc4 -np5'
	extras='-na0.0 -P/xs -nc4 -np5'
c	extras='-na0.0 -nc1 -np5 -hmax500 -hmin490 -P/null'

	do i=1,7
	  rhomax=3.*sqrt(pi*(hmax)/2./period(i)/light)*180./pi
	  be_max=rhomax+3.*sqrt(rhomax)
	  be_max=rhomax
	  do j=1,ngen(i)
c choose alpha from 10 to 170
 5	    alpha=ran1(iseed)*180.
            if(alpha.lt.10. .or. alpha.gt.170.)goto 5
            alpha=45.
c choose beta from -be_max to +be_max
 10	    beta=(2*ran1(iseed)-1.0)*be_max
	    if(abs(beta).gt.be_max)goto 10

	    i_alpha=nint(alpha)
	    i_beta=nint(beta)
            write(10,'('' rm -f fort.66'')')
	    if(i_beta.lt.0)then
	      write(10,'('' beam '',a,'' -pd'',f4.2,'' -al'',i3.3,
     +          '' -be-'',i2.2,'' -r'',i10.10)')
     +          extras,period(i),i_alpha,abs(i_beta),iseed
	    else
	      write(10,'('' beam '',a,'' -pd'',f4.2,'' -al'',i3.3,
     +          '' -be'',i2.2,'' -r'',i10.10)')
     +          extras,period(i),i_alpha,i_beta,iseed
	    endif 
            write(10,'('' cat fort.66 >> res '')')
	  enddo
	  write(*,*)i,period(i),be_max
	enddo

	close(10)

	end
