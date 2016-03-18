      subroutine profparams(noise,totflux)

      implicit none

      include 'useful.inc'

      integer i, zpos
      real hmax_low,hmax_high
      real totflux,totlin,totcirc
      real fwhm10,noise,snr
      character angle*10,ach*20,bch*20,kch*20

      profymin=1.e6
      profymax=0.0
      totflux=0.0
      totlin=0.0
      totcirc=0.0
      zpos=0

c compute the total and linear flux
      do i=1,np
	 totflux=totflux+profi(i)
         if (profi(i).gt.profymax) then
           profymax=profi(i)
           zpos=i
         endif
         if (profi(i).lt.profymin) profymin=profi(i)
         if (profv(i).lt.profymin) profymin=profv(i)
	 totlin=totlin+profl(i)
	 totcirc=totcirc+profv(i)
c         write(67,*) xprof(i),profi(i),profl(i)
      enddo
      totflux=totflux/real(np)
      totlin=totlin/real(np)
      totcirc=totcirc/real(np)
      snr=totflux*sqrt(float(np))/noise
      if(profymin.gt.-3.0*noise)profymin=-3.0*noise

c compute the FWHM
      hmax_low=0.		
      hmax_high=0.		
      do i = zpos,1,-1
	if(profi(i) .gt. profymax/2.)hmax_low=xprof(i)
      enddo
      do i = zpos,np
	if(profi(i) .gt. profymax/2.)hmax_high=xprof(i)
      enddo
      fwhm=hmax_high-hmax_low

c	write(*,*)profymax,zpos,xprof(zpos),fwhm,hmax_high,hmax_low
c compute the FWHM10
      hmax_low=0.		
      hmax_high=0.		
      do i = zpos,1,-1
	if(profi(i) .gt. profymax/10.)hmax_low=xprof(i)
      enddo
      do i = zpos,np
	if(profi(i) .gt. profymax/10.)hmax_high=xprof(i)
      enddo
      fwhm10=hmax_high-hmax_low

      if(verb)then
        write(*,*)
        write(*,'('' Max width   = '',f10.2)')widthmax
        write(*,'('' Total flux  = '',f10.2)')totflux
        write(*,'('' Peak  flux  = '',f10.2)')profymax
        write(*,'('' FWHM        = '',f10.2)')fwhm
        write(*,'('' FW 10%      = '',f10.2)')fwhm10
        write(*,'('' SNR         = '',f10.2)')snr
        write(*,'('' Linear frac = '',f10.2)')totlin/totflux*100.
        write(*,'('' Circ frac   = '',f10.2)')totcirc/totflux*100.
        if(key(1:1) .eq. '-')then
	  write(*,'('' NOT DETECTED !!! '')')
        else
          write(*,'('' Class       = '',a1)')key
          if(fwhm .gt. widthmax)then
            write(*,'('' FWHM is more than max width !!! '')')
          endif
        endif
      endif
      
      return
      end

