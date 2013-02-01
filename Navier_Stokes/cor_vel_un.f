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

      subroutine cor_vel_un(us, u_no_pres, An,
     &                      dpdx, density,
     &                      nx, ny, dx, dy, dt)

c calculates u^** and v^** i.e. eq, 12-13
c k is a flag: 1 for u, 1 for v
c see pg. 13 for notes



      implicit double precision (a-h, o-z)

      double precision us(nx,ny) ! u_star or v_star
      double precision u_no_pres(nx,ny) 

      double precision An(nx,ny) ! A^* times ones

      double precision dpdx(nx,ny) ! gradiend of pressure ps

      a = (dx*dy)/dt    

      do j = 1,ny
        do i = 1,nx

          us(i,j) = u_no_pres(i,j)
     &      - dx * dy * dpdx(i,j) / (density * (a - 0.5d0 * An(i,j) ))

        enddo
      enddo

      return
      end

 
