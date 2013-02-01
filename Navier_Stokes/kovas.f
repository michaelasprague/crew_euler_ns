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

      subroutine kovas(u, v, p, xmin, ymin, dx, dy,
     &                 nx, ny, t, viscosity, density)

! calculates Taylor Green solution at time t


      implicit double precision (a-h, o-z)


      double precision u(nx,ny), v(nx,ny), p(nx,ny)
 
      pi = acos(-1.d0)
      Re = 1.d0/viscosity
      ! lambda
      c = 0.5d0*Re - ( 0.25d0*Re**2 + 4.d0*pi**2 )**0.5d0


      umax = 0.d0
      vmax = 0.d0

      do j = 1,ny

        do i = 1, nx

          x = xmin+0.5d0*dx +(i-1)*dx
          y = ymin+0.5d0*dy +(j-1)*dy

          u(i,j) = 1.d0 - exp(c*x) * cos(2.d0*pi*y)

          v(i,j) =  c * exp(c*x) * sin(2.d0*pi*y) / (2.d0 * pi)
      
          if (abs(u(i,j)) .gt. umax) umax = abs(u(i,j))
          if (abs(v(i,j)) .gt. vmax) vmax = abs(v(i,j))
   
          p(i,j) = ( 1.d0 - exp(2.d0*c* x)) / 2.d0

          write(79,*) x, p(i,j)
                         
        enddo

      enddo

c     write(*,*) 'umax,vmax',umax,vmax
c     stop

      return
      end

    
