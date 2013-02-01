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

      subroutine max_speeds_euler(u, v, p, r, u_max, sound_max, nx, ny)

      implicit double precision(a-h,o-z)

! solution at time t 
      double precision u (nx, ny)
      double precision v (nx, ny)
      double precision r (nx, ny)
      double precision p (nx, ny)

      write(*,*) 'nx = ', nx
      write(*,*) 'ny = ', ny
      write(*,*) 'u(1,1) = ', u(1,1)
      write(*,*) 'v(1,1) = ', v(1,1)
      write(*,*) 'p(1,1) = ', p(1,1)
      write(*,*) 'r(1,1) = ', r(1,1)

      sound_max = 0.d0
      u_max = 0.d0
      do j = 1, ny
        do i = 1, nx
          u_mag = sqrt(u(i,j)**2 + v(i,j)**2)
          if (u_mag .gt. u_max) u_max = u_mag
          if (p(i,j) .lt. 0.) stop 'p negative'
          if (r(i,j) .lt. 0.) stop 'r negative'
          sound = sqrt(p(i,j)/r(i,j))
          if (sound .gt. sound_max) sound_max = sound
        enddo
      enddo

      return
      end
