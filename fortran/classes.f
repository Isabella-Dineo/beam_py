      subroutine classify(np,compmax,xprof,profi,derivi,noise,verb,
     +                    singles,doubles,multis,nons,key)

      implicit none
      integer np,compmax,i,n_max

      integer loc(compmax),nons,singles,doubles,multis,dip
      real xprof(np),profi(np),derivi(np)
      real cmaxi,cdip,noise
      character key*1
      logical verb

      dip=0
      cmaxi=0.0
      cdip=0.0
c find the maxima in the profile
      n_max=0
      do i=1,np-2
c       write(*,*)i,xprof(i),cmaxi,cdip,profi(i),derivi(i)
        if (profi(i).ge.3.0*noise) then 
          if (n_max.eq.0) then
            if (derivi(i).ge.0.0.and.derivi(i+1).lt.0.0.and.
     $          derivi(i+2).lt.0.0)then
c             if (derivi(i).gt.0.0.and.derivi(i+1).lt.0.0)then
                n_max=n_max+1
                loc(n_max)=xprof(i)
                cmaxi=profi(i)
                dip=0
              end if
            else
              if (derivi(i).ge.0.0.and.derivi(i+1).lt.0.0.and.
     $            derivi(i+2).lt.0.0)then
                if (dip.eq.1.and.profi(i).gt.cdip+1.0*noise) then
                  n_max=n_max+1
                  loc(n_max)=xprof(i)
                  cmaxi=profi(i)
                  dip=0
                end if
              else
                 if (derivi(i).lt.0.0.and.derivi(i+1).ge.0.0.and
     $                .derivi(i+2).ge.0) then
                  if (profi(i).lt.cmaxi-1.0*noise.and.dip.eq.0) then
                    dip=1
                    cdip=profi(i)
                  end if
                end if
              end if
            end if
          else
             if (derivi(i).lt.0.0.and.derivi(i+1).ge.0.0.and.derivi(i+2)
     $            .ge.0) then
              if (profi(i).lt.cmaxi-1.0*noise.and.dip.eq.0) then
                dip=1
                cdip=profi(i)
              end if
            end if
          end if
               
      end do

c      if(verb)then
c        do i=1,n_max
c          write(*,*)'Maxima at:', loc(i)
c        end do
c      endif
      key="-"
      if (n_max.eq.0) then
        key='-'
        nons=nons+1
        write(*,*)'nons are ',nons
      end if
      if (n_max.eq.1) then
        singles=singles+1
        key='S'
        write(*,*)'singles are ',singles
      end if
      if (n_max.eq.2) then
        doubles=doubles+1
        key='D'
        write(*,*)'doubles are ',doubles
      end if
      if (n_max.gt.2) then
        multis=multis+1
        key='M'
        write(*,*)'multis are ',multis
      end if
            
      return
      end
