      subroutine plotlos(xlos,ylos,np)

      implicit none
      integer i,j,np
      real xlos(np),ylos(np)

      call pgsci(2)
      call pgslw(3)
      call pgmove(xlos(1),ylos(1))
      do i=2,np
        call pgdraw(xlos(i),ylos(i))
      enddo

c draw evenly spaced dots every 20 degrees
      do i=1,18
	j=nint(float(i)/18.*float(np))
	call pgpt(1,xlos(j),ylos(j),17)
      enddo

      call pgsci(1)
      call pgslw(1)

      end

