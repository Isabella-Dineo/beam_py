C *************************************************************************
      program fits2epn

c     Reads a PARKES FITS file and writes an EPN file
c     The header is read mostly following the sequence
c     of the epn-header variables (not optimal :-)


      implicit none
      include 'epnhdr.inc'

      integer status,unit,readwrite,blocksize,hdutype,ntable
      integer felem,nelems,nullj,diameter,nfound,irow,colnum
      integer nrows,ncols,icol
      integer iargc,epnrecs
      integer ibl,ipl,blkn
      integer j
      integer*2 tbi,m3da(1024,1025,4)
      integer*4 fbi,ifr


      real fbr
      real*8 dr,ofs,sint
      real*8 dar(1025)
      real*8 chbw,realstart

      real nulle,density
      character ifile*40,nullstr*1,name*8,ttype(21)*10
      character decsign
      character*2 polc(4)
      character istring80*80,istring24*24,istring16*16,istring8*8
      character*8 xname,vrn,chu
      character*80 oname

      logical anynull
      logical there


c     Read and respond to Arguments
      
      if (iargc().lt.1) then
         write(*,*)''
         write(*,*)'fits2epn --- (last update 12/2/2003)'
         write(*,*)'by Aris Karastergiou'
         write(*,*)'email: aris@physics.usyd.edu.au'
         write(*,*)''
         write(*,*)'USAGE:  fits2epn <FITS-file>'
         write(*,*)''
         write(*,*)'Produces .epn file'
         write(*,*)''
         stop
      endif
      
      call getarg(1,ifile)
      
      inquire(file=ifile,exist=there)
      if(.not.there) then
         write(*,*)''
         write(*,*)'Cannot find file ',ifile
         stop
      endif


C     The STATUS and other parameters must always be initialized.

      oname='test.epn'
      status=0
      irow=1
      felem=1
      nelems=1
      nullstr=' '
      nullj=0
      nulle=0.
      
C     Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)
      
C     Open the FITS file
      readwrite=0
      call ftopen(unit,ifile,readwrite,blocksize,status)
      
      
C     Find out main things
c     How many subints are in the file: HDU SUBINT, rows...

      xname='SUBINT'
      hdutype=-1
      call ftmnhd(unit,hdutype,xname,0,status)
      call ftgnrw(unit,nrows,status)
      write(*,*)'Fits file has ',nrows,' sub-integrations.'
      epnrecs=nrows

c     Main header details
      history='Parkes FITS file turned EPN'

c     Pulsar name
      xname='PSREPHEM'
      call ftmnhd(unit,hdutype,xname,0,status)
      vrn='PSR_NAME'
      call FTGCNO(unit,1,vrn,icol,status)
      call ftgcvs(unit,icol,irow,felem,nelems,nullstr,istring80,anynull
     $     ,status)
      write(*,*)istring80
      write(jname,'(''J'',a11)')istring80
      write(*,*)'Your source is PSR ',jname
      cname=jname
c***********************************************************************
c     CAREFUL period dm rm are taken from PSREPHEM!!! MUST CHANGE
c***********************************************************************
      xname='PSREPHEM'
      call ftmnhd(unit,hdutype,xname,0,status)

c     Pbar
      vrn='IF0'
      call FTGCNO(unit,1,vrn,icol,status)
      call ftgcvj(unit,icol,irow,felem,nelems,nullstr,fbi,anynull
     $     ,status)
      vrn='FF0'
      call FTGCNO(unit,1,vrn,icol,status)
      call ftgcvd(unit,icol,irow,felem,nelems,nullstr,dr,anynull
     $     ,status)
      pbar=1./(fbi+dr)
      write(*,*)'Period:',pbar

c     DM
      vrn='DM'
      call FTGCNO(unit,1,vrn,icol,status)
      call ftgcvd(unit,icol,irow,felem,nelems,nullstr,dr,anynull
     $     ,status)
      dm=dr
      write(*,*)'DM:',dm

c     RM
      vrn='RM'
      call FTGCNO(unit,1,vrn,icol,status)
      call ftgcvd(unit,icol,irow,felem,nelems,nullstr,dr,anynull
     $     ,status)
      rm=dr
      write(*,*)'RM:',rm

c     Other obscure EPN variables

      catref=''
      bibref=''

c     RA and Dec

      vrn='RAJ'
      call FTGCNO(unit,1,vrn,icol,status)
      call ftgcvs(unit,icol,irow,felem,nelems,nullstr,istring24,anynull
     $     ,status)
      read(istring24(1:2),'(i2)')rah
      read(istring24(4:5),'(i2)')ram
      read(istring24(7:12),'(f6.0)')ras
      write(*,*)rah,ram,ras

      vrn='DECJ'
      call FTGCNO(unit,1,vrn,icol,status)
      call ftgcvs(unit,icol,irow,felem,nelems,nullstr,istring24,anynull
     $     ,status)
      decsign=istring24(1:1)
      read(istring24(2:3),'(i2)')ded
      read(istring24(5:6),'(i2)')dem
      read(istring24(8:13),'(f6.0)')des
      if(decsign.eq.'-')ded=-ded
      write(*,*)ded,dem,des
      
c     Telescope name
      ntable=1 !the name is guaranteed to be in the first header
      call ftmahd(unit,1,hdutype,status)
      call ftgkys(unit,'TELESCOP',istring24,istring80,status)
      read(istring24(1:8),'(a8)')telname
      write (*,*) 'Observed in ',telname
c     Telescope x,y,z
      call ftgkyd(unit,'ANT_X   ',dr,istring80,status)
      xtel=dr
      call ftgkyd(unit,'ANT_Y   ',dr,istring80,status)
      ytel=dr
      call ftgkyd(unit,'ANT_Z   ',dr,istring80,status)
      ztel=dr
      write(*,*)'Which is located at x,y,z:',xtel,ytel,ztel

c     Date of observation

      call ftgkys(unit,'STT_DATE',istring24,istring80,status)
      read(istring24(1:4),'(i4)')cdy
      read(istring24(6:7),'(i2)')cdm
      read(istring24(9:10),'(i4)')cdd
      write (*,*)'On ',cdd,'/',cdm,'/',cdy

c     Epoch
      call ftgkyj(unit,'STT_IMJD',fbi,istring80,status)
      epoch=real(fbi)
      write(*,*)'MJD:',epoch

c     Get start time of first subint while here...
      call ftgkyj(unit,'STT_SMJD',fbi,istring80,status)
      realstart=real(fbi)
      write(*,*)realstart,' seconds after the start of the epoch.'



c     OPOS
      xname='SUBINT'
      hdutype=-1
      call ftmnhd(unit,hdutype,xname,0,status)
      vrn='POS_ANG'
      call FTGCNO(unit,1,vrn,icol,status)
      call ftgcve(unit,icol,irow,felem,nelems,nullstr,fbr,anynull
     $     ,status)
      opos=fbr
      write(*,*)'OPOS:',opos

c     PAFLAG and  TIMFLAG
      paflag=' '
      timflag='U'

c     scanno and subscan

      scanno=1
      subscan=1

c     Number of polarizations
      xname='HISTORY'
      call ftmnhd(unit,hdutype,xname,0,status)
      vrn='NPOL'
      call FTGCNO(unit,1,vrn,icol,status)
      call ftgcvi(unit,icol,irow,felem,nelems,nullstr,tbi,anynull
     $     ,status)
      npol=tbi
      write(*,*)'Your file contains ',npol,' polarizations,'

c     Polarization identifier

      vrn='POL_TYPE'
      call FTGCNO(unit,1,vrn,icol,status)
      call ftgcvs(unit,icol,irow,felem,nelems,nullstr,istring80,anynull
     $     ,status)
      read(istring80(1:2),'(a2)')polc(1)
      read(istring80(3:4),'(a2)')polc(2)
      read(istring80(5:6),'(a2)')polc(3)
      read(istring80(7:8),'(a2)')polc(4)
      write(*,*)'of this type: ',polc(1),polc(2),polc(3),polc(4)


c     Frequency bands per polarization

      vrn='NCHAN'
      call FTGCNO(unit,1,vrn,icol,status)
      call ftgcvi(unit,icol,irow,felem,nelems,nullstr,tbi,anynull
     $     ,status)
      nfreq=tbi
      write(*,*)'each containing ',nfreq,' frequency channels,'

c     and, while you're at it, get the channel bandwidth
c     need to get units as well. CAREFUL!!!!

      vrn='CHAN_BW'
      call FTGCNO(unit,1,vrn,icol,status)
      call ftgcvd(unit,icol,irow,felem,nelems,nullstr,dr,anynull
     $     ,status)
      chbw=abs(dr)
      write(*,*)'of a bandwidth of ',chbw,' MHz each.'
      

c     Number of phase bins per period

      vrn='NBIN_PRD'
      call FTGCNO(unit,1,vrn,icol,status)
      call ftgcvi(unit,icol,irow,felem,nelems,nullstr,tbi,anynull
     $     ,status)
      nbin=tbi
      write(*,*)'Each period is sampled into ',nbin,' bins.'

c     Sampling interval

      vrn='TBIN'
      call FTGCNO(unit,1,vrn,icol,status)
      call ftgcvd(unit,icol,irow,felem,nelems,nullstr,dr,anynull
     $     ,status)
      tbin=dr*1000000.
      write(*,*)'Time per bin:',tbin,' usec.'

c     Number of integrated pulses per data block!CAREFUL this changes
c     from one SUBINT to the next!!!!

      xname='SUBINT'
      call ftmnhd(unit,hdutype,xname,0,status)
      vrn='TSUBINT'
      call FTGCNO(unit,1,vrn,icol,status)
      call ftgcvd(unit,icol,irow,felem,nelems,nullstr,dr,anynull
     $     ,status)
      nint=int(dr/pbar)
      write(*,*)nint

c     A bit about the cal needs to go in here

      ncal=1
      lcal=0
      
c     TRES = TBIN

      tres=tbin

c     Fluxflag is still unset

      fluxflag=' '

c     End of main EPN header

c***********************************************************************
c     Start of each subheader. Loop around the polarizations
c     and the frequency channels...

      sint=0.0                  !Time offset of first block

      do ibl=1,epnrecs
         blkn=0

c     Read in MONSTER 3d data array!!
         xname='SUBINT'
         call ftmnhd(unit,hdutype,xname,0,status)
         vrn='DATA'
         call FTGCNO(unit,1,vrn,icol,status)
         irow=ibl
         nelems=4198400
         call ftgcvi(unit,icol,irow,felem,nelems,nullstr,m3da,anynull
     $        ,status)
         write(*,*)'Read in data...'

c     Read in center frequencies for the freq-channels

         call ftmnhd(unit,hdutype,xname,0,status)
         vrn='DAT_FREQ'
         call FTGCNO(unit,1,vrn,icol,status)
         irow=ibl
         nelems=nfreq           !Total array elements to read
         call ftgcvd(unit,icol,irow,felem,nelems,nullstr,dar
     $        ,anynull,status)
         nelems=1               !Set it back to normal 1
         
         
         do ipl=1,npol !Loop around the polarizations
            do ifr=1,nfreq !Loop around the frequency channels
               blkn=blkn+1      !current sub-block number
               
c     The identification TAG: pol tag 
               write(idfield(blkn),'(a2)')polc(ipl)
               
c     The ordinal frequency channel number for each polarization
               nband(blkn)=ifr
               
c     Number of channels averaged into current channel
               navg(blkn)=1
               
c     Centre sky frequency of current channel

               f0(blkn)=dar(ifr)/1000.
               f0u(blkn)='GHz'
               
c     Channel Bandwidth (has been assigned a value further up)
               
               df(blkn)=chbw
               dfu(blkn)='MHz'
               
c     Time stamp of current subint

               tstart(blkn)=(realstart+sint)

c     Parallactic angle info !FOR NOW HARDWIRED

               plflag(blkn)='A'
               plval(blkn)=0.0

c     OTHER EPN STUFF

               scale(blkn)=0.0
               offset(blkn)=0.0
               rms(blkn)=0.0
               papp(blkn)=pbar !SOME CARE MUST BE TAKEN HERE

c     The dreaded DATA!!!

               do j=1,nbin
                  rawdata(blkn,j)=real(m3da(j,ifr,ipl))
               end do
            end do
         end do
         call rwepn(oname,+2,ibl,.false.)
         write(*,*)'Time stamp of Subint ',ibl
         write(*,'(f17.5)')100000*tstart(1)
c      Increment time stamp!
         xname='SUBINT'
         call ftmnhd(unit,hdutype,xname,0,status)
         vrn='OFFS_SUB'
         call FTGCNO(unit,1,vrn,icol,status)
         irow=ibl
         call ftgcvd(unit,icol,irow,felem,nelems,nullstr,dr
     $        ,anynull,status)
         ofs=dr
         vrn='TSUBINT'         
         call FTGCNO(unit,1,vrn,icol,status)
         irow=ibl
         call ftgcvd(unit,icol,irow,felem,nelems,nullstr,dr
     $        ,anynull,status)
         sint=sint+dr
      end do

      end
