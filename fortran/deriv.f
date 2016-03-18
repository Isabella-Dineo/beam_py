      subroutine deriv(dim,inar,outar)

c     subroutine to calculate the central difference
c     as an approximation of the derivative
c
c     arguments:
c
c     inar     : the input array
c     outar    : the output array
c     dim      : the dimension of the two arrays
c
c     outar(1) and outar(dim) will use forward and backward difference
c     respectively

      integer i,j,dim,lag
      real inar(dim),outar(dim)
      real jo

c      outar(1)=inar(2)-inar(1)

      do i=3,dim-2
         outar(i)=0.25*(inar(i+2)-inar(i-2))
      end do
      outar(2)=0.5*(inar(3)-inar(1))
      outar(1)=outar(2)
      outar(dim-1)=0.5*(inar(dim)-inar(dim-2))
      outar(dim) = inar(dim)-inar(dim-1)
      end
         
