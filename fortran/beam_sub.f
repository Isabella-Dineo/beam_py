      subroutine beam_sub(al,be,hmin,hmax,ncomp,npatch,cw,noise, period
     $     ,freq,iseed,same_again,flatpol,dobeam,doprof,doab,dm,flux)

      implicit none

      include 'useful.inc'

      real*8 al, be, phi, rhomax, rhomax2
      real*8 wtmp, pwmax2
      real*8 x1, x2, dx
      real*8 y1, y2, dy
      real*8 si(nz,nz),sq(nz,nz),su(nz,nz),sv(nz,nz)

      integer i, j, k, iseed,iseed1
      integer ncomp,npatch,xind,yind,npatch_core
      integer total,detected
      real derivi(np),pa_rvm(np),tmp_pa
      real ran1,period,noise,rlc,freq,dm
      real peaki,sn,flux,aveIbin
      real*8 hmin,hmax,hrange
      real*8 cw,rho_core,phi0,psi0,al_rvm,be_rvm
      real*8 ab_time,ab_deg,ab_xoff(compmax),ab_yoff(compmax)
      real*8 ab_xoff3(compmax),ab_yoff3(compmax)
      logical same_again,dobeam,doprof,flatpol,doab
      real*8 sinD, asinD
      character xopt*10,yopt*10

c Initialise parameters
      x1=-180.
      x2=180.
      dx=(x2-x1)/float(nz)
      y1=-180.
      y2=180.
      dy=(y2-y1)/float(nz)
c setup the longitude array
      phi=-180.
      do i=1,np
         phi=phi+360./float(np)
         xprof(i)=phi
      enddo
c     
      do i=1,nz
        do j=1,nz
          si(i,j)=0.0
          sq(i,j)=0.0
          su(i,j)=0.0
          sv(i,j)=0.0
        enddo
      enddo
      
c     light cylinder radius and maximum rho and height range
      rlc=period*light/2./pi
      rhomax=3.*sqrt(pi*(hmax)/2./period/light)*180./pi
      hrange=hmax-hmin   
    
      if (verb) then 
         write(*,'('' Pulsar period  = '',f10.6)')period
         write(*,'('' Light cylinder = '',f10.1)')rlc
         write(*,'('' Min height     = '',f10.1)')hmin
         write(*,'('' Max height     = '',f10.1)')hmax
         write(*,'('' Height range   = '',f10.1)')hrange
         write(*,'('' Number of components = '',i4)')ncomp
         write(*,'('' Number of patches    = '',i4)')npatch
         write(*,'('' Observing frequency  = '',f6.2)')freq
      end if 
      
c     compute the component heights
      if(.not.same_again)then
         do i=1,ncomp
           hcomp(i)=hmin+hrange*(ran1(iseed)**2)
         enddo
      endif
      
      do i=1,ncomp
         
c     compute the component heights (freq dependent) and hence rho
           new_hcomp(i)=0.6*hcomp(i)*((freq/1.4)**-0.667)
     +                     + 0.4*hcomp(i)
           rho(i)=3.*sqrt(pi*new_hcomp(i)/2./period/light)*180./pi

c aberration stuff - aberration time just given by h_em/c
c convert to degrees and compute the change in the location
c of the pole on the x,y beam plane
	   if(doab)then
	     ab_time = new_hcomp(i)/light
             ab_deg = ab_time/period*360.
             call mapphi(al,0.0,ab_deg,ab_xoff(i),ab_yoff(i))
             call mapphi(al,0.0,3.*ab_deg,ab_xoff3(i),ab_yoff3(i))
c       write(*,*)ab_deg,ab_xoff(i),ab_yoff(i),ab_xoff3(i),ab_yoff3(i)
           else
	     ab_xoff(i)=0.0
	     ab_yoff(i)=0.0
	     ab_xoff3(i)=0.0
	     ab_yoff3(i)=0.0
           endif

c     decide what to do with pw !!
c      pw(i)=1. + rho(i)/5.
c     pw(i)=rho(i)/5.+0.5
c     pw(i)=sqrt(rho(i))
c     Mitra and Rankin Style + Aris jitter
c           pw(i)=0.08*2.45*sqrt(new_hcomp(i)/10.0/period)+1.
           pw(i)=0.2*2.45*sqrt(new_hcomp(i)/10.0/period)

      enddo
            
      rhomax2=0.0
      pwmax2=0.0
      do i=1,ncomp
        if (verb) then
           write(*,'('' Component '',i2)')i
           write(*,'('' Height of component is '',f12.6)'
     $                 )new_hcomp(i)
           write(*,'('' Rho   for component is '',f12.6)')rho(i)
           write(*,'('' Width for this component is '',f12.6)')
     +            pw(i)
        end if
	if(rho(i).gt.rhomax2)then
	  rhomax2=rho(i)
	  pwmax2=pw(i)
	endif
        call patchbeam(rho(i),npatch,pw(i),iseed,x1,dx,y1,dy
     $             ,si,sq,su,sv,nz,pcd,same_again,flatpol,
     $              ab_xoff(i),ab_yoff(i),ab_xoff3(i),ab_yoff3(i),
     $              new_hcomp(i))
       enddo

c     add the core if necessary
         if(cw.gt.0.)then
           rho_core=0.0
           npatch_core=1
	   call patchbeam(rho_core,npatch_core,cw,iseed,x1,dx,y1,dy,si
     $           ,sq,su,nz,pcd,same_again,flatpol)
         endif
            
c     compute the maximum profile width
c     see equation 3.26 in Lorimer and Kramer
c     added in extra pw/2 to account for grazing
         wtmp=sinD(rhomax2/2.+pwmax2/2.)**2-sinD(be/2.)**2
         if(wtmp.lt.0. .or. al+be.eq.0.)then
            if(verb)write(*,'('' Warning: beta > rho !! '')')
              widthmax=360.
            else
              widthmax=sqrt(wtmp/sinD(al)/abs(sinD(al+be)))
              if(widthmax.gt.1)then
                widthmax=360.
              else
                widthmax=4.*asin(widthmax)*180./pi
              endif
         endif

         profxmax=widthmax/2.+20.
         if(profxmax.gt.180.)profxmax=180.
         profxmin=-profxmax
         if(verb)write(*,'('' Max prof width = '',f10.1)')widthmax
               
c     compute the line of sight
c     this returns the x,y coords as a function of observer longitude
c     and there is thus a direct mapping between observer longitude
c     and the particular magnetic field line
         call getlos(al,be,xlos,ylos,thetalos,np)
               
c get the rvm curve - pass angles in degrees
c this is the angle between the los and the projected field line
         phi0=0.0
         psi0=0.0
	 al_rvm=al
	 be_rvm=be
         call rvm(al_rvm,be_rvm,phi0,psi0,xprof,pa_rvm,np)
               
c     create the profiles in I, Q, U from the beam and line of sight
         do j=1,np
           xind=nint((xlos(j)-x1)/dx)
           yind=nint((ylos(j)-y1)/dy)
           k=np-j+1
           profi(k)=si(xind,yind)
           profq(k)=sq(xind,yind)
           profu(k)=su(xind,yind)
           profl(k)=sqrt(profq(k)**2+profu(k)**2)
           profv(k)=sv(xind,yind)
         enddo
cc **************************************************
cc This is where a call to scattering routine should be made
cc **************************************************
         call scatter(profi,np,dm,period,freq,verb)
         call scatter(profq,np,dm,period,freq,verb)
         call scatter(profu,np,dm,period,freq,verb)         

         call scatter(profv,np,dm,period,freq,verb)
c     find the maxima in the profile and classify
         call deriv(np,profi,derivi)
         call classify(np,compmax,xprof,profi,derivi,noise,verb,key)
c     normalise for S/N
         peaki=0.0
         do i=1,np
            if (profi(i).gt.peaki) peaki = profi(i)
         end do
         if (noise.gt.0)then
            sn=noise
            noise=1.0
         else
            sn=1
         end if
         do i=1,np
            profi(i)=profi(i)*sn/peaki
            profq(i)=profq(i)*sn/peaki
            profu(i)=profu(i)*sn/peaki
            profv(i)=profv(i)*sn/peaki
         end do
c get profile info before adding noise
	 call profparams(noise,aveIbin) ! aveIbin is the mean I per bin

c     add noise if required
         iseed1=-int(secnds(0.0))
         if(noise.gt.0.)then
            call addnoise(iseed1,profi,np,noise)
            call addnoise(iseed1,profq,np,noise)
            call addnoise(iseed1,profu,np,noise)
            call addnoise(iseed1,profv,np,noise)
         endif
c     Now the stokes profiles have noise of 1 and peak flux of S/N
c     rescale to get the flux right:
         flux=flux/aveIbin
         peaki=peaki*flux
         profymax=profymax*flux
         profymin=profymin*flux
         noise=noise*flux
         do i=1,np
            profi(i)=profi(i)*flux
            profq(i)=profq(i)*flux
            profu(i)=profu(i)*flux
            profv(i)=profv(i)*flux
         end do
c     compute the final PA swing from the Q,U profiles
c     need to debias L with the noise
c     note the mapping between the actual PA of the field line (thetalos)
c     and the RVM curve (pa_rvm)
         do i=1,np-1
            profl(i)=sqrt(profq(i)**2+profu(i)**2) - sqrt(2.)*noise
            tmp_pa=0.0
            epa(i)=0.0
            if (profl(i).gt.3.*noise)then
              tmp_pa=0.5*atan2(profu(i),profq(i))*180./pi
	      epa(i)=noise/profl(i)*180./pi
            endif
	    if(tmp_pa.eq.0.0)then
              profpa(i)=0.0
            else
              profpa(i)=tmp_pa-thetalos(i)+pa_rvm(i)
c	      profpa(i)=180./pi*acos(((xlos(i+1)-xlos(i))*xlos(i)
c     +          +(ylos(i+1)-ylos(i))*ylos(i))/(sqrt((xlos(i+1)
c     +          -xlos(i))**2+(ylos(i+1)-ylos(i))**2)*sqrt(xlos(i)**2
c     +          +ylos(i)**2)))-90.0
 20	      if(profpa(i).lt.-90.)then
                profpa(i)=profpa(i)+180.
                goto 20
              endif
 21	      if(profpa(i).gt.90.)then
                profpa(i)=profpa(i)-180.
                goto 21
              endif
            endif
c           write(55,*)xprof(i),tmp_pa,thetalos(i),profpa(i),pa_rvm(i)
         enddo

c     plot things, including the beam, the profile and the PA swing
         xopt='bcnst'
         yopt='bcnst'
	 if(dobeam)then
           call plotbeam(si,nz,profxmin,profxmax)
           call plotlos(xlos,ylos,np)
	 endif
	 if(doprof)then
           call pgsvp(0.57,0.99,0.10,0.70)
           call plotprof(al,be,xopt,yopt)
           call plotpa(pa_rvm)
	 endif
               
      end
