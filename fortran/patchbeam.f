c     subroutine to fill a beam at particular rho with patches of
c     emission, and output I(x,y) Q(x,y) and U(x,y) where x,y are the
c     coordinates on the plane

      subroutine patchbeam(rho,npatch,pw,iseed,x1,dx,y1,dy,si,sq,su,sv
     $     ,nz,pcd,samemag,flatpol,ab_xoff,ab_yoff,ab_xoff3,ab_yoff3,
     $     new_hcomp)

c     rho :   rho of this ring (input)
c     npatch: number of patches
c     pw  :   patch width in sigma
c     iseed:  ran1 seed
c     si  :   I(x,y) (output)
c     sq  :   Q(x,y) (output)
c     su  :   U(x,y) (output)
c     nz  :   array dimension (input)
c     x1,y1,dx,dy: map array indices to x and y


      implicit none
      integer nz,i,j,k,npatch,iseed
      real*8 si(nz,nz), sq(nz,nz), su(nz,nz), sv(nz,nz)
      real*8 rho,x1,dx,y1,dy,x,y,pw
      real*8 pcx(npatch),pcy(npatch), pcd(npatch), pi
      real*8 iamp(npatch),angle,frp
      real*8 radius1,result,frpol(npatch),xangle,yangle,result2
      real*8 ia,qa,ua,xtmp,ytmp,frcirc(npatch)
      real*8 ab_xoff,ab_yoff,ab_xoff3,ab_yoff3,new_hcomp
      real ran1
      parameter (pi=3.1415926545)
      logical samemag,flatpol


c     find the centres of all the patches in x, y coordinates
      if(npatch.gt.0)then
      do i=1, npatch
         if (.not. samemag) pcd(i) = 2.*pi*ran1(iseed) ! patch centre in rad
         pcx(i)=rho*sin(pcd(i))
         pcy(i)=rho*cos(pcd(i))
         iamp(i)=10.0

c compute fractional polarization at random
c have a chance for an orthogonal mode flip
c         frpol(i)=iamp(i)*ran1(iseed)
c         frpol(i)=0.1+new_hcomp/500.0
          frpol(i)=1.0
c         if(ran1(iseed).lt.abs(0.8-frpol(i))) frpol(i)=-frpol(i)
         if(ran1(iseed).lt.0.5) frpol(i)=-frpol(i)
	 frpol(i)=frpol(i)*iamp(i)
c compute fractional circular polarization 
c between -25 and +25 %
	  frcirc(i)=iamp(i)*(ran1(iseed)-0.5)/2.

c         write(*,*) pcd(i)*180.0/pi, pcx(i),pcy(i),frpol(i)
      end do

      do i=1,nz
         x=x1+float(i)*dx
         do j=1,nz
            y=y1+float(j)*dy
	    angle=atan2(x+ab_xoff/3.,y+ab_yoff/3.)
            do k=1,npatch

c total intensity with aberration offset
	       xtmp=(x-pcx(k)-ab_xoff)*(x-pcx(k)-ab_xoff)
	       ytmp=(y-pcy(k)-ab_yoff)*(y-pcy(k)-ab_yoff)
               radius1=sqrt(xtmp+ytmp)
	       if(radius1/pw .gt. 3.)then
	         result=0.0
                 result2=0.0
	       else
                 result=exp(-(((radius1)/pw)**2))
                 result2=exp(-(((radius1)/(0.8*pw))**2))
c the I flux depends how far you are away from the patch up to 3 sigma
                 si(i,j)=si(i,j)+iamp(k)*result
                 sv(i,j)=sv(i,j)+frcirc(k)*result
c now compute Q and U given the location of the patch
c the position angle corresponds to the particular field line
c at the centre of the patch
	         if(flatpol)then
                   sq(i,j)=sq(i,j)+frpol(k)*result2*cos(2.*pcd(k))
                   su(i,j)=su(i,j)+frpol(k)*result2*sin(2.*pcd(k))
	         else
c or compute Q and U with a position angle
c corresponding to the actual field line
                   sq(i,j)=sq(i,j)+frpol(k)*result2*cos(2.*angle)
                   su(i,j)=su(i,j)+frpol(k)*result2*sin(2.*angle)
	         endif
               endif
            end do
         end do
       end do

c if npatch=0 make a uniformly illuminated cone
       else
         ia=10.0
         frp=ia*ran1(iseed)
         frp=ia
         if(ran1(iseed).lt.0.0)frp=-frp
	 do k=1,360,int(pw)+1
           angle=real(k)*pi/180.
c offset the angles by the aberration values
           xangle=rho*sin(angle)+ab_xoff
           yangle=rho*cos(angle)+ab_yoff
           do j=1,nz
             y=y1+float(j)*dy
	     ytmp=(y-yangle)*(y-yangle)
             do i=1,nz
               x=x1+float(i)*dx
               xtmp=(x-xangle)*(x-xangle)
               radius1=sqrt(xtmp+ytmp)
	       if(radius1/pw .gt. 3.)then
	         result=0.0
	       else
                 result=exp(-(((radius1)/pw)**2))
                 si(i,j)=si(i,j)+ia*result
c offset the angles by the aberration values
	         angle=atan2(x+ab_xoff3,y+ab_yoff3)
                 sq(i,j)=sq(i,j)+frp*result*cos(2.*angle)
                 su(i,j)=su(i,j)+frp*result*sin(2.*angle)
	       endif
	     enddo
	   enddo
	 enddo
       endif


      return
      end

