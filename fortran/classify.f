      subroutine classify(np,compmax,xprof,profi,derivi,noise,verb,key)

      implicit none
      integer np,compmax,i,n_max

      integer loc(compmax),dip
      real xprof(np),profi(np),derivi(np)
      real cmaxi,cdip,noise
      character key*10
      logical verb

      dip=0
      cmaxi=0.0
      cdip=0.0
c find the maxima in the profile
      n_max=0
      do i=1,np-2
c       write(*,*)i,xprof(i),cmaxi,cdip,profi(i),derivi(i)
        if (profi(i).ge.5.0*noise) then 
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
      end if
      if (n_max.eq.1) then
        key='S'
      end if
      if (n_max.eq.2) then
        key='D'
      end if
      if (n_max.gt.2) then
        key='M'
      end if
            
      return
      end
