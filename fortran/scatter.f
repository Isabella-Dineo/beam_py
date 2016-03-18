      subroutine scatter(prof,np,dm,period,freq,verb)

      implicit none
      integer i,tb,tbf,np
      real prof(np),tprof(np),dm,period,freq
      real conv(np)
      real*8 tbin,tau
      real tcrit, bincrit,scr
      logical verb

      if (dm.eq.0) goto 20
c     scr is the scatter response
c     conv is the convolution result
      do i=1,np
         conv(i)=0.0
      end do

c     first compute the time of each bin in ms
      tbin=period*1000/np
c     now find the scattering time given the dm and frequency
      tau=-6.46+0.154*log10(dm)+1.07*(log10(dm))**2-3.86*log10(freq)
      tau=10**tau
      if (verb) write(*,*)'Scattering:' ,tau, dm,period,freq,tbin
c     try to calculate what sort of delays are meaningful for tau
c     0.01=exp(-tcrit/tau) <=> log(0.01)/-tcrit=1/tau
      tcrit=-tau/log(0.01)
      bincrit=tcrit/tbin
      if (bincrit.gt.10*np) bincrit=10*np
      if (bincrit.lt.np) bincrit=np
      if (verb) write(*,*)'bincrit:',bincrit
c     Now follow the recipe for convolution:
c     1. Transpose the profile:
      do i=1,np
         tprof(np-i+1)=prof(i)
      end do
c     2. Slide along the time axis for as long as the scatter response
c     is still high
      tb=1                      ! straight bin number
      tbf=1                     ! folded bin number
      do while (tb.lt.bincrit+np)
         do i=1,np
            if (tb+i-2.gt.np) then
               scr=exp(-(tb+i-2-np)*tbin/tau)
            else
               scr=0
            end if
            conv(tbf)=conv(tbf)+tprof(i)*scr
         end do
         tb=tb+1
         tbf=tbf+1
         if (tbf.gt.np) tbf=1
      end do
      do i=1,np
         prof(i)=conv(i)
      end do
 20    continue
      end

               
            
