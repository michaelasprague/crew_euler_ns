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

      subroutine divergence(div, u, v,
     &                      u_bc_ew, 
     &                      v_bc_ns,
     &                      div_rms, div_max,
     &                      dx, dy, nx, ny)

      implicit double precision(a-h,o-z)

      double precision div(nx,ny), u(nx,ny), v(nx,ny)
      double precision u_bc_ew(ny,2)
      double precision v_bc_ns(nx,2)

      do j = 2, ny-1
        do i = 2, nx-1
          div(i,j) = 0.5*dy*(u(i+1,j) - u(i-1,j))
     &             + 0.5*dx*(v(i,j+1) - v(i,j-1))
        enddo
      enddo

      do i = 2, nx-1
        j = 1 !south
        div(i,j) = 0.5*dy*(u(i+1,j) - u(i-1,j))
     &           + 0.5*dx*(v(i,j+1) + v(i,j)) 
     &           - dx * v_bc_ns(i,1)

        j = ny ! north
        div(i,j) = 0.5*dy*(u(i+1,j) - u(i-1,j))
     &           + dx * v_bc_ns(i,2)
     &           - 0.5 * dx * (v(i,j) + v(i,j-1)) 

      enddo

      do j = 2, ny-1
          i = 1 ! west
          div(i,j) = 
     &           + 0.5*dy*(u(i+1,j) + u(i,j))
     &           - dy*u_bc_ew(j,1)
     &           + 0.5*dx*(v(i,j+1) - v(i,j-1))

          i = nx ! east
          div(i,j) = 
     &           + dy*u_bc_ew(j,2)
     &           - 0.5*dy*(u(i,j) + u(i-1,j))
     &           + 0.5*dx*(v(i,j+1) - v(i,j-1))
      enddo

      i = 1  ! west
      j = 1  ! south
      div(i,j) =
     &           + 0.5*dy*(u(i+1,j) + u(i,j))
     &           - dy*u_bc_ew(j,1)
     &           + 0.5*dx*(v(i,j+1) + v(i,j)) 
     &           - dx*v_bc_ns(i,1)

      i = nx ! east
      j = 1  ! south
      div(i,j) =
     &           + dy*u_bc_ew(j,2)
     &           - 0.5*dy*(u(i,j) + u(i-1,j))
     &           + 0.5*dx*(v(i,j+1) + v(i,j)) 
     &           - dx*v_bc_ns(i,1)

      i = 1  ! west
      j = ny ! north
      div(i,j) =
     &           + 0.5*dy*(u(i+1,j) + u(i,j))
     &           - dy*u_bc_ew(j,1)
     &           + dx*v_bc_ns(i,2)
     &           - 0.5*dx*(v(i,j) + v(i,j-1)) 

      i = nx ! east
      j = ny ! north
      div(i,j) =
     &           + dy*u_bc_ew(j,2)
     &           - 0.5*dy*(u(i,j) + u(i-1,j))
     &           + dx*v_bc_ns(i,2)
     &           - 0.5*dx*(v(i,j) + v(i,j-1)) 


      div_max = 0.d0
      div_rms = 0.d0
      do j = 1, ny
        do i = 1, nx
          div(i,j) = abs(div(i,j))
          if (div(i,j) .gt. div_max) div_max = div(i,j)
          div_rms = div_rms + div(i,j)**2
        enddo
      enddo

      div_rms = sqrt(div_rms / float(nx*ny))

      return
      end
