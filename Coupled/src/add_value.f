      subroutine add_value(u,val, nx,ny)
 

c adds specified scalar value to given array u

      implicit double precision (a-h,o-z)
      double precision u (nx, ny)

      do i = 1, nx 
        do j = 1, ny

          u(i,j) = u(i,j) + val

        enddo
      enddo

      return
      end
