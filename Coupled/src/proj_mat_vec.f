!***********************************************************************
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory (NREL)
!
!    This is free software: you can redistribute it and/or modify it
!    under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the
!    License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License 
!    along with this program.
!    If not, see <http://www.gnu.org/licenses/>.
!
!***********************************************************************
!    This code was created at NREL by Michael A. Sprague and Ignas 
!    Satkauskas  and was meant for open-source distribution.
!
!    Software was created under funding from a Shared Research Grant 
!    from the Center for Research and Education in Wind (CREW), during
!    the period 01 October 2011 - 31 January 2013.
!
!    http://crew.colorado.edu/ 
!
!    Questions?  Please contact Michael Sprague:
!    email:  michael.a.sprague@nrel.gov
!
!***********************************************************************

      subroutine proj_mat_vec(w, p, nx, ny, dx, dy)

      implicit double precision (a-h,o-z)

      double precision p(nx,ny), w(nx,ny)


      dx2 = dx*dx
      dy2 = dy*dy

c Inside

      do j = 2, ny-1
        do i = 2, nx-1
           w(i,j) = 
     &            + (p(i+1,j) - p(i,j)  )/dx2
     &            - (p(i,j)   - p(i-1,j))/dx2
     &            + (p(i,j+1) - p(i,j)  )/dy2
     &            - (p(i,j)   - p(i,j-1))/dy2

        enddo
      enddo


c South

      j = 1
      do i = 2, nx-1

!       dpdy_s = -(2.*p(i,j) - 3.*p(i,j+1) + p(i,j+2))/dy
!       dpdy_n = (-p(i,j) + p(i,j+1))/dy

        dpdy_s = (2.*p(i,j) )/dy
        dpdy_n = (-p(i,j) + p(i,j+1))/dy

        w(i,j) = 
     &            + (p(i+1,j) - p(i,j)  )/dx2
     &            - (p(i,j)   - p(i-1,j))/dx2
     &            + dpdy_n / dy
     &            - dpdy_s / dy

      enddo

c North

      j = ny 
      do i = 2, nx-1

!       dpdy_s = (+p(i,j) - p(i,j-1))/dy
!       dpdy_n = -(-2.*p(i,j) + 3.*p(i,j-1) - p(i,j-2))/dy

        dpdy_s = (+p(i,j) - p(i,j-1))/dy
        dpdy_n = (-2.*p(i,j) )/dy

        w(i,j) = 
     &            + (p(i+1,j) - p(i,j)  )/dx2
     &            - (p(i,j)   - p(i-1,j))/dx2
     &            + dpdy_n / dy
     &            - dpdy_s / dy

      enddo

c West

      i = 1
      do j = 2, ny-1

!       dpdx_w = -(2.*p(i,j) - 3.*p(i+1,j) + p(i+2,j))/dx
!       dpdx_e = (-p(i,j) + p(i+1,j))/dx
        dpdx_w = (2.*p(i,j) )/dx
        dpdx_e = (-p(i,j) + p(i+1,j))/dx

        w(i,j) = 
     &            + dpdx_e / dx
     &            - dpdx_w / dx
     &            + (p(i,j+1) - p(i,j)  )/dy2
     &            - (p(i,j)   - p(i,j-1))/dy2

      enddo

c East

      i = nx
      do j = 2, ny-1

!       dpdx_e = -(-2.*p(i,j) + 3.*p(i-1,j) - p(i-2,j))/dx
!       dpdx_w = (p(i,j) - p(i-1,j))/dx
        dpdx_e = (-2.*p(i,j) )/dx
        dpdx_w = (p(i,j) - p(i-1,j))/dx

        w(i,j) = 
     &            + dpdx_e / dx
     &            - dpdx_w / dx
     &            + (p(i,j+1) - p(i,j)  )/dy2
     &            - (p(i,j)   - p(i,j-1))/dy2

      enddo


c Corners:

c SW
      i = 1
      j = 1
!     dpdx_w = -(2.*p(i,j) - 3.*p(i+1,j) + p(i+2,j))/dx
!     dpdy_s = -(2.*p(i,j) - 3.*p(i,j+1) + p(i,j+2))/dy
      dpdx_w = (2.*p(i,j) )/dx
      dpdy_s = (2.*p(i,j) )/dy
      w(i,j) =  0.d0
     &            + (p(i+1,j) - p(i,j)  )/dx2
     &            - dpdx_w / dx
     &            + (p(i,j+1) - p(i,j)  )/dy2
     &            - dpdy_s / dy

c SE
      i = nx 
      j = 1
!     dpdx_e = -(-2.*p(i,j) + 3.*p(i-1,j) - p(i-2,j))/dy
!     dpdy_s = -(2.*p(i,j) - 3.*p(i,j+1) + p(i,j+2))/dy
      dpdx_e = (-2.*p(i,j) )/dx
      dpdy_s = (2.*p(i,j) )/dy
      w(i,j) =  0.d0
     &            + dpdx_e / dx
     &            - (p(i,j)   - p(i-1,j))/dx2
     &            + (p(i,j+1) - p(i,j)  )/dy2
     &            - dpdy_s / dy

c NW
      i = 1
      j = ny
!     dpdy_n = -(-2.*p(i,j) + 3.*p(i,j-1) - p(i,j-2))/dy
!     dpdx_w = -(2.*p(i,j) - 3.*p(i+1,j) + p(i+2,j))/dx
      dpdy_n = (-2.*p(i,j) )/dy
      dpdx_w = (2.*p(i,j) )/dx
       w(i,j) =  0.d0
     &            + (p(i+1,j) - p(i,j)  )/dx2
     &            - dpdx_w / dx
     &            + dpdy_n /dy
     &            - (p(i,j)   - p(i,j-1))/dy2

c NE
      i = nx
      j = ny
!     dpdy_n = -(-2.*p(i,j) + 3.*p(i,j-1) - p(i,j-2))/dy
!     dpdx_e = -(-2.*p(i,j) + 3.*p(i-1,j) - p(i-2,j))/dx
      dpdy_n = (-2.*p(i,j) )/dy
      dpdx_e = (-2.*p(i,j) )/dx
      w(i,j) =  0.d0
     &            + dpdx_e / dx
     &            - (p(i,j)   - p(i-1,j))/dx2
     &            + dpdy_n / dy
     &            - (p(i,j)   - p(i,j-1))/dy2


      return
      end

