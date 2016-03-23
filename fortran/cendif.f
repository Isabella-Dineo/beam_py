      subroutine cendif(inar,outar,dim,jo,lag)

c     subroutine to calculate the central difference
c     as an approximation of the derivative
c
c     arguments:
c
c     inar     : the input array
c     outar    : the output array
c     dim      : the dimension of the two arrays
c     jo       : takes in the desired threshold 
c     lag      : writes out the first lag 
c
c     outar(1) and outar(dim) will use forward and backward difference
c     respectively

      integer i,j,dim,lag
      real inar(dim),outar(dim)
      real jo

c      outar(1)=inar(2)-inar(1)

      do i=2,dim-1
         outar(i)=0.5*(inar(i+1)-inar(i-1))
      end do
      
      outar(1)=outar(2)
      outar(dim) = inar(dim)-inar(dim-1)

      do i=2,dim-1
         if ((outar(i).gt.jo.and.outar(i+1).gt.jo)) then
            lag=i
            goto 10
         end if
      end do
 10   continue
      end
         
