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

      subroutine get_kovas_pbc(u_bc_ew, u_bc_ns, v_bc_ew, v_bc_ns,
     &                        dp_bc_ew, dp_bc_ns, 
     &                        xmin, xmax, ymin, ymax,
     &                        t, nx, ny, dx, dy, viscosity)

      implicit double precision (a-h,o-z)

      double precision u_bc_ew(ny,2), u_bc_ns(nx,2)
      double precision v_bc_ew(ny,2), v_bc_ns(nx,2)
      double precision dp_bc_ew(ny,2), dp_bc_ns(nx,2)

      ! bc_ew( number of elements, face 1-west,2-east, velocity component)
      ! bc_ns (number of elements, face 1-south, 2 north)

      pi = acos(-1.d0)
      Re = 1.d0 / viscosity
      c = 0.5d0*Re - ( 0.25d0*Re**2.d0 + 4.d0*pi**2.d0 )**0.5d0

c     write(*,*) 'viscosity = ', viscosity
c     write(*,*) 'pi = ', pi
c     write(*,*) 'Re = ', Re
c     write(*,*) 'c (lambda) = ', c

      if (abs(c + 0.963741) .gt. 0.0001) stop 'problem in get_kovas_pbc'


!  u(x,y) = 1 - E**(c*x)*Cos(2*Pi*y)
!  v(x,y) = (c*E**(c*x)*Sin(2*Pi*y))/(2.*Pi)

! south   (  -1/2 < x < 1, y = -1/2 )
! north  (  -1/2 < x < 1, y = 1+1/2 )
      do i = 1, nx

        x = xmin + 0.5d0*dx +float(i-1)*dx

        y = ymin
        u_bc_ns(i,1)  = 1.d0 - exp(c*x) * cos(2.d0 * pi * y)
        v_bc_ns(i,1)  = c * exp(c*x) * sin(2.d0*pi*y)/ (2.d0 * pi)
        dp_bc_ns(i,1) = 0.d0 

        y = ymax
        u_bc_ns(i,2)  = 1.d0 - exp(c*x) * cos(2.d0 * pi * y)
        v_bc_ns(i,2)  = c * exp(c*x) * sin(2.d0*pi*y)/ (2.d0 * pi)
        dp_bc_ns(i,2) = 0.d0 

      enddo


! west  ( x = -1/2, -1/2 < y < 1+1/2 )
! east  ( x = 1,    -1/2 < y < 1+1/2 )
      do j = 1, ny

        y = ymin + 0.5d0*dy + float(j-1)*dy

        x = xmin
        u_bc_ew(j,1)  = 1.d0 - exp(c*x) * cos(2.d0*pi*y)
        v_bc_ew(j,1)  = c * exp(c*x) * sin(2.d0*pi*y) / (2.d0 * pi)
        dp_bc_ew(j,1) = -c*exp(2.d0*c*x)

        x = xmax
        u_bc_ew(j,2)  = 1.d0 - exp(c*x) * cos(2.d0*pi*y)
        v_bc_ew(j,2)  = c * exp(c*x) * sin(2.d0*pi*y)/ (2.d0 * pi)
        dp_bc_ew(j,2) = -c * exp(2.d0*c*x)

      enddo

      return
      end
