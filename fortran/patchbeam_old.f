c     subroutine to fill a beam at particular rho with patches of
c     emission, and output I(x,y) Q(x,y) and U(x,y) where x,y are the
c     coordinates on the plane

      subroutine patchbeam(rho, npatch, pw, iseed, x1,dx,y1,dy,si,sq,su
     $     ,nz)

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
      real*8 rho, si(nz,nz), sq(nz,nz), su(nz,nz),pw
      real*8 x1,dx, y1,dy,x,y
      real*8 pcx(npatch),pcy(npatch), pcd(npatch), pi
      real*8 patchi,patchq,patchu
      real ran1
      parameter (pi=3.1415)

c     find the centres of all the patches in x, y coordinates

      do i=1, npatch
         pcd(i) = 2.*pi*ran1(iseed) ! patch centre in rad
         pcx(i)=rho*sin(pcd(i))
         pcy(i)=rho*cos(pcd(i))
c         write(*,*) pcd(i)*180.0/pi, pcx(i),pcy(i)
      end do

      do i=1,nz
         x=x1+float(i)*dx
         do j=1,nz
            y=y1+float(j)*dy
            do k=1,npatch
               si(i,j)=si(i,j)+patchi(x,y,pw,pcx(k),pcy(k))
               sq(i,j)=sq(i,j)+patchq(x,y,pw,pcx(k),pcy(k),pcd(k))
               su(i,j)=su(i,j)+patchu(x,y,pw,pcx(k),pcy(k),pcd(k))
            end do
         end do
      end do

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      real*8 function patchi(x,y,pw,pcx,pcy)
      implicit none
      real*8 rho,x,y,r,pw,pcx,pcy
      real*8 result

      r=sqrt((x-pcx)*(x-pcx)+(y-pcy)*(y-pcy))

      result=10.0*exp(-(((r)/pw)**2))
      patchi=result
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      real*8 function patchq(x,y,pw,pcx,pcy,pcd)
      implicit none
      real*8 rho,x,y,r,pw,q,pcx,pcy,pcd
      real*8 result

      q=10.0/sqrt(2.)*cos(pcd)

      r=sqrt((x-pcx)*(x-pcx)+(y-pcy)*(y-pcy))

      result=q*exp(-(((r)/pw)**2))
      patchq=result
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      real*8 function patchu(x,y,pw,pcx,pcy,pcd)
      implicit none
      real*8 rho,x,y,r,pw,u,pcx,pcy,pcd
      real*8 result

      u=10.0/sqrt(2.)*sin(pcd)

      r=sqrt((x-pcx)*(x-pcx)+(y-pcy)*(y-pcy))
      result=u*exp(-(((r)/pw)**2))
      patchu=result
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
