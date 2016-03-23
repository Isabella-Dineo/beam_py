      program beam2fits
      
      implicit none
      include 'useful.inc'

      integer i,j,nargs,imjd,mjds
      real rdum,period,freq
      real*8 ddum,mjdf
      character beamfile*20, fitstemp*20,fitsfile*24
      character cdum*24
      nargs = IARGC()

      if (nargs.lt.8) then
         write(*,*)'Usage: beam2fits <beamfile> <fitstemp> <fitsfile> 
     &  <period(s)> <frequency(GHz)> <the 3 MJDs>'
         stop
      end if
      call getarg(1,beamfile)
      call getarg(2,fitstemp)
      call getarg(3,fitsfile)
      call getarg(4,cdum)
      read(cdum(1:),*)period
      call getarg(5,cdum)
      read(cdum(1:),*)freq
c The MJD stuff
      call getarg(6,cdum)
      read(cdum(1:),*)imjd
      call getarg(7,cdum)
      read(cdum(1:),*)mjds
      call getarg(8,cdum)
      read(cdum(1:),*)mjdf

c     Read in beamfile and fill up arrays
      open (unit=77, file=beamfile,form='formatted', status='unknown')
      do i=1,np
         read(77,*)profi(i),profq(i),profu(i),profv(i)
      end do
      close(77)
      call beamfits(period,freq,fitstemp,fitsfile,imjd,mjds,mjdf)
      end
      
