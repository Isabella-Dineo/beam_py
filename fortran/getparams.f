	subroutine getparams(al,be,hmin,hmax,ncomp,npatch,cw,noise,
     1   period,freq,plotdev,iseed,verb,flatpol,dobeam,doprof,doab,dm,
     2   outfile,flux)

	implicit none
	integer i,ncomp,npatch,nargs,iseed
	real*8 al,be,hmin,hmax,width,cw
	real noise,period,freq,dm,flux
	character string*50,plotdev*50,outfile*20
	logical happy,verb,dobeam,doprof,flatpol,doab

	nargs = IARGC()

c defaults
	happy=.true.
	al=45.
	be=5.
	iseed=0
	hmax=1000.
	hmin=20.
	ncomp=4
	npatch=10
	cw=0.0
	noise=0.0
	flux=1.0
        period=0.1
	freq=1.4
	verb=.false.
        flatpol=.false.
	dobeam=.true.
	doprof=.true.
        doab=.false.
	plotdev='/xs'
	outfile='simIQUV.dat'
	dm=0
	if(nargs .eq. 0)then
	  write(*,'('' Required options '')')
	  write(*,'('' -r<iseed> '')')
	  write(*,*)
	  write(*,'('' Optional options '')')
	  write(*,'('' -al<alpha> (def='',f6.2,'')'')')al
	  write(*,'('' -be<beta> (def='',f6.2,'')'')')be
	  write(*,'('' -P<plot device> (def=/xs)'')')
	  write(*,'('' -A<ascii file> (def=simIQUV.dat)'')')
	  write(*,'('' -hmax<max height> (def='',f8.2,'')'')')hmax
	  write(*,'('' -hmin<max height> (def='',f6.2,'')'')')hmin
	  write(*,'('' -nc<number of components> (def='',i4,'')'')')ncomp
	  write(*,'('' -np<number of patches> (def='',i4,'')'')')npatch
c	  write(*,'('' -pw<patch width> (def='',f6.2,'')'')')pw
	  write(*,'('' -co<core width> (def=no core)'')')
	  write(*,'('' -na<peak S/N> (def=no noise)'')')
	  write(*,'('' -S<continuum eq. flux> (def=peak flux 1)'')')
	  write(*,'('' -pd<pulse period> (def=0.1)'')')
	  write(*,'('' -fr<frequency in GHz> (def='',f6.2,'')'')')freq
	  write(*,'('' -dm<the DM> (def='',f6.2,'')'')')dm
	  write(*,'('' -flat<flat polarization> (false)'')')
	  write(*,'('' -ab<do aberration> (false)'')')
	  write(*,'('' -v (verbose)'')')
	  write(*,*)
	endif

	do i=1,nargs
	  call getarg(i,string)

c check for seed option 
	  if (index(string,'-r') .gt. 0)then
	     read(string(3:),*)iseed
	     if(iseed.gt.0)iseed=-iseed
	     call ran1(iseed)
	  endif

c check for alpha option 
	  if (index(string,'-al') .gt. 0)then
	     read(string(4:),*)al
	  endif

c check for beta option 
	  if (index(string,'-be') .gt. 0)then
	     read(string(4:),*)be
	  endif

c check for hmax option 
	  if (index(string,'-hmax') .gt. 0)then
	     read(string(6:),*)hmax
	  endif

c check for hmin option 
	  if (index(string,'-hmin') .gt. 0)then
	     read(string(6:),*)hmin
	  endif

c check for number of components option 
	  if (index(string,'-nc') .gt. 0)then
	     read(string(4:),*)ncomp
	     if(ncomp.lt.1)ncomp=1
	  endif

c check for number of patches option 
	  if (index(string,'-np') .gt. 0)then
	     read(string(4:),*)npatch
	  endif

c check for patch width option 
c	  if (index(string,'-pw') .gt. 0)then
c	     read(string(4:),*)pw
c	  endif

c check for core option 
	  if (index(string,'-co') .gt. 0)then
	     read(string(4:),*)cw
	  endif

c check for noise option 
	  if (index(string,'-na') .gt. 0)then
	     read(string(4:),*)noise
	  endif

c check for flux option 
	  if (index(string,'-S') .gt. 0)then
	     read(string(3:),*)flux
	  endif

c check for period option 
	  if (index(string,'-pd') .gt. 0)then
	     read(string(4:),*)period
	  endif

c check for frequency option 
	  if (index(string,'-fr') .gt. 0)then
	     read(string(4:),*)freq
	  endif

c check for flatpol
	  if (index(string,'-flat') .gt. 0)then
	     flatpol=.true.
	  endif

c check for aberration
	  if (index(string,'-ab') .gt. 0)then
	     doab=.true.
	  endif
	  
c check for DM option 
	  if (index(string,'-dm') .gt. 0)then
	     read(string(4:),*)dm
	  endif
	  
c check for verbosity
	  if (index(string,'-v') .gt. 0)then
	     verb=.true.
	  endif
	  
c check for plot device option 
	  if (index(string,'-P') .gt. 0)then
	     read(string(3:),'(a)')plotdev
	     if(plotdev.eq.'/null' .or. plotdev .eq.'/NULL')then
               dobeam=.false.
               doprof=.false.
	     endif
	  endif
c check for ascii output option
	  if (index(string,'-A') .gt. 0)then
	     read(string(3:),'(a)')outfile
	  endif
	  
	enddo

	if(iseed.eq.0)then
	  write(*,'('' You must specify a value for the seed !! '')')
	  happy=.false.
	endif

        if(.not. happy)then
          write(*,*)
	  write(*,'('' Rerun with appropriate parameters !!'')')
	  stop
	endif

	return
	end

