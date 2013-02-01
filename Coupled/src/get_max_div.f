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

      subroutine get_max_div(div_max,bd_flux,u,v,dx,dy,nx,ny)

      implicit double precision (a-h,o-z)

      double precision u (nx, ny-1)
      double precision v (nx-1, ny)
 
      div_max = 0.d0

      do j = 1, ny-1

        do i = 1, nx-1

          div = dy*( u(i+1,j) - u(i,j)) + dx * ( v(i,j+1) - v(i,j))

          if (abs(div) .gt. div_max) div_max = abs(div)

        enddo

      enddo


      flux = 0.d0
! south
      j = 1
      do i = 1,nx-1
        flux = flux + v(i,j)*dx
      enddo
! north
      j = ny
      do i = 1,nx-1
        flux = flux - v(i,j)*dx
      enddo
! west
      i = 1 
      do j = 1,ny-1
        flux = flux + u(i,j)*dy
      enddo
! east
      i = nx
      do j = 1,ny-1
        flux = flux - u(i,j)*dy
      enddo

      bd_flux = flux
!     write(*,*) 'boundary flux = ',flux

      return
      end
