      subroutine beamfits(period,freq,template,outfile,imjd,mjds,mjdf)

      implicit none

      include 'useful.inc'

      integer frch
      parameter (frch=1)
      integer i,j
      integer unit, status,idum,cn
      integer imjd,mjds
      real*8 tbin,f0,coefs(5),ddum,bw_ch,dfreq,mjdf
      real freq,period,iquv(np,frch,4),bw,mjd
      real dum(frch),dum2(frch,4)
      character template*14,date*24,cdum*24,comment*24
      character outfile*24

c     Observing frequency of sim data:
      freq=1000.0*freq
c     Bandwidth
      bw=256.0
      bw_ch=256.0/frch
c     Open the fits file to write
      status=0
c      template='testtemp.txt'
      call ftgiou(unit,status)
      call fttplt(unit, outfile, template,status)
      write(*,*)'outfile and template opened'
c     go to primary HDU
      call ftmahd(unit,1,idum,status)
      call ftgstm(date,idum,status)
      call ftgkey(unit,"DATE",cdum,comment ,status)
      call ftukys(unit,"DATE",date,comment,status)
c      call ftukye(unit,"BMAJ",1.0,9,"",status)
c      call ftukye(unit,"BMIN",1.0,9,"",status)
c      call ftukye(unit,"SCANLEN",1.0,9,"",status)
c      call ftukye(unit,"STT_OFFS",0.0,9,"",status)
c      call ftukye(unit,"STT_LST",0.0,9,"",status)
c      call ftukye(unit,"OBSFREQ",freq,9,"",status)
c      call ftukye(unit,"OBSBW",bw,4,"",status)
c      call ftukyj(unit,"OBSNCHAN",frch,"",status)
      call ftukyj(unit,"STT_IMJD",imjd,"",status)
      call ftukyj(unit,"STT_SMJD",mjds,"",status)
      call ftukyd(unit,"STT_OFFS",mjdf,9,"",status)
c-----------------------------------------
c-----------------------------------------
c     and so on for the header
c     Now go to the subint stuff
      call ftmnhd(unit,2,"SUBINT",0,status)
      tbin=period/np
      call ftukyd(unit,"TBIN",tbin,9,"",status)
      call ftukyj(unit,"NBIN",np,"",status)
      call ftukyj(unit,"NCHAN",frch,"",status)
      call ftukyd(unit,"CHAN_BW",bw_ch,3,"",status)
c     Change the dimensions of the various arrays
      ddum=60.0
      call ftpcld(unit,1,1,1,1,ddum,status) !tsubint
      ddum=0.0
      call ftpcld(unit,2,1,1,1,ddum,status) !subint offset
      ddum=5000.0
      call ftpcld(unit,3,1,1,1,ddum,status) !dummy lst
      ddum=0.0
      call ftpcld(unit,4,1,1,1,ddum,status) !dummy RA
      call ftpcld(unit,5,1,1,1,ddum,status) !dummy dec
      call ftpcle(unit,8,1,1,1,0.0,status)! feed angle
      call ftpcle(unit,9,1,1,1,0.0,status)! pos angle
      call ftpcle(unit,10,1,1,1,0.0,status)! para angle
      call ftpcle(unit,11,1,1,1,0.0,status)! telescope az
      call ftpcle(unit,12,1,1,1,0.0,status)! telescope zen
      call ftmvec(unit,13,frch,status)! data freqs
      call ftmvec(unit,14,frch,status)! data weights
      call ftmvec(unit,15,4*frch,status)! data offsets
      call ftmvec(unit,16,4*frch,status)! data scales
      do j=1,frch
         dum(j)=freq-bw/2+(j-1)*bw_ch
      end do
      call ftpcle(unit,13,1,1,frch,dum,status)
      do j=1,frch
         dum(j)=1.0
      end do
      call ftpcle(unit,14,1,1,frch,dum,status)
      do i=1,4
         do j=1,frch
            dum2(j,i)=0.0
         end do
      end do
      call ftpcle(unit,15,1,1,4*frch,dum2,status)
      do i=1,4
         do j=1,frch
            dum2(j,i)=1.0
         end do
      end do
      call ftpcle(unit,16,1,1,4*frch,dum2,status)
c     and so on
c 
c     now change the dimensions of the data array
       idum=np*4*frch! nbins*npol*nchan
c      idum=np! nbins*npol*nchan
      call ftmvec(unit,17,idum,status)
      call ftukys(unit,"TDIM17","(512,1,4)","",status)
c      call ftukys(unit,"TDIM17","(512)","",status)
      do j=1,frch
         do i=1,np
            iquv(i,j,1)=profi(i)
         end do
         do i=1,np
            iquv(i,j,2)=profq(i)
         end do
         do i=1,np
            iquv(i,j,3)=profu(i)
         end do
         do i=1,np
            iquv(i,j,4)=profv(i)
         end do
      end do
c      idum=idum*4
      call ftpcle(unit,17,1,1,idum,iquv,status)
c     call ftpcle(unit,17,1,1,np,profi,status)
c-----------------------------------
c-----------------------------------
c     Go to the polyco stuff
      do i=1,5
         coefs(i)=0.0
      end do
      mjd=real(imjd)+real(mjds)/86400
      call ftmnhd(unit,2,"POLYCO",0,status)
      call ftgcno(unit,.true.,"REF_F0",cn,status)
      f0=1.0/period
      call ftpcld(unit,cn,1,1,1,f0,status)
      call ftgcno(unit,.true.,"REF_MJD",cn,status)
      call ftpcle(unit,cn,1,1,1,mjd,status)
      call ftgcno(unit,.true.,"REF_PHS",cn,status)
      call ftpcle(unit,cn,1,1,1,0.0,status)
      call ftgcno(unit,.true.,"PRED_PHS",cn,status)
      call ftpcle(unit,cn,1,1,1,0.0,status)
      call ftgcno(unit,.true.,"LGFITERR",cn,status)
      call ftpcle(unit,cn,1,1,1,-6.0,status)
      call ftgcno(unit,.true.,"COEFF",cn,status)
      call ftpcle(unit,cn,1,1,5,coefs,status)
      call ftgcno(unit,.true.,"NCOEF",cn,status)
      call ftpclj(unit,cn,1,1,1,5,status)
      call ftgcno(unit,.true.,"NSPAN",cn,status)
      call ftpclj(unit,cn,1,1,1,1440,status)
      call ftgcno(unit,.true.,"NPBLK",cn,status)
      call ftpclj(unit,cn,1,1,1,1,status)
      call ftgcno(unit,.true.,"REF_FREQ",cn,status)
      call ftpcle(unit,cn,1,1,1,1420.0,status)
      call ftgcno(unit,.true.,"NSITE",cn,status)
      call ftpcls(unit,cn,1,1,1,"0",status)
      call ftgcno(unit,.true.,"DATE_PRO",cn,status)
      call ftpcls(unit,cn,1,1,1,date,status)
      call ftgcno(unit,.true.,"POLYVER",cn,status)
      call ftpcls(unit,cn,1,1,1,"Splc",status)
c-----------------------------------
c-----------------------------------
c     Go to the history stuff
      call ftmnhd(unit,2,"HISTORY",0,status)
      call ftgcno(unit,.true.,"DATE_PRO",cn,status)
      call ftpcls(unit,cn,1,1,1,date,status)
      call ftgcno(unit,.true.,"PROC_CMD",cn,status)
      call ftpcls(unit,cn,1,1,1,"Simulated data",status)
      write(*,*)'-----------'
      call ftgcno(unit,.true.,"SCALE",cn,status)
      call ftpcls(unit,cn,1,1,1,"Jy",status)
      call ftgcno(unit,.true.,"POL_TYPE",cn,status)
      call ftpcls(unit,cn,1,1,1,"IQUV",status)
      call ftgcno(unit,.true.,"NSUB",cn,status)
      call ftpclj(unit,cn,1,1,1,1,status)
      call ftgcno(unit,.true.,"NPOL",cn,status)
      call ftpclj(unit,cn,1,1,1,4,status)
      call ftgcno(unit,.true.,"NBIN",cn,status)
      call ftpclj(unit,cn,1,1,1,np,status)
      call ftgcno(unit,.true.,"NBIN_PRD",cn,status)
      call ftpclj(unit,cn,1,1,1,np,status)
      call ftgcno(unit,.true.,"TBIN",cn,status)
      call ftpcld(unit,cn,1,1,1,tbin,status)
      call ftgcno(unit,.true.,"CTR_FREQ",cn,status)
      dfreq=freq
      call ftpcld(unit,cn,1,1,1,dfreq,status)
      call ftgcno(unit,.true.,"NCHAN",cn,status)
      call ftpclj(unit,cn,1,1,1,frch,status)
      call ftgcno(unit,.true.,"CHAN_BW",cn,status)
      call ftpcld(unit,cn,1,1,1,bw_ch,status)
      call ftgcno(unit,.true.,"NBIN",cn,status)
      call ftpclj(unit,cn,1,1,1,np,status)
      call ftgcno(unit,.true.,"PAR_CORR",cn,status)
      call ftpclj(unit,cn,1,1,1,0,status)
      call ftgcno(unit,.true.,"FA_CORR",cn,status)
      call ftpclj(unit,cn,1,1,1,0,status)
      call ftgcno(unit,.true.,"RM_CORR",cn,status)
      call ftpclj(unit,cn,1,1,1,0,status)
      call ftgcno(unit,.true.,"DEDISP",cn,status)
      call ftpclj(unit,cn,1,1,1,0,status)
      call ftgcno(unit,.true.,"DDS_MTHD",cn,status)
      call ftpcls(unit,cn,1,1,1,"NONE",status)
      call ftgcno(unit,.true.,"SC_MTHD",cn,status)
      call ftpcls(unit,cn,1,1,1,"off",status)
      call ftgcno(unit,.true.,"CAL_MTHD",cn,status)
      call ftpcls(unit,cn,1,1,1,"NONE",status)
      call ftgcno(unit,.true.,"CAL_FILE",cn,status)
      call ftpcls(unit,cn,1,1,1,"NONE",status)
      call ftgcno(unit,.true.,"RFI_MTHD",cn,status)
      call ftpcls(unit,cn,1,1,1,"A",status)
      call ftgcno(unit,.true.,"IFR_MTHD",cn,status)
      call ftpcls(unit,cn,1,1,1,"NONE",status)
      call ftclos(unit, status)
      call ftfiou(unit, status)      
      end
