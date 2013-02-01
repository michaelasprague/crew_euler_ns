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

      subroutine inflow_bcs(u_bc_ew, u_bc_ns, v_bc_ew, v_bc_ns,
     &                        dp_bc_ew, dp_bc_ns, 
     &                         nx, ny, dx, dy, viscosity,flux_tot)

c inflow bc conditions

      implicit double precision (a-h,o-z)

      double precision u_bc_ew(ny,2), u_bc_ns(nx,2)
      double precision v_bc_ew(ny,2), v_bc_ns(nx,2)
      double precision dp_bc_ew(ny,2), dp_bc_ns(nx,2)

      ! bc_ew( number of elements, face 1-west,2-east, velocity component)
      ! bc_ns (number of elements, face 1-south, 2 north)

      !vel=flux_tot/(2.d0*nx*dx+2.d0*ny*dy)

      vel = 1.d0

      write(*,*) 'vel',vel

      do j = 1, ny
       
       u_bc_ew(j,1) = vel !west
       v_bc_ew(j,1) = vel !0.d0
c       dp_bc_ew(j,1) = 0.d0

       u_bc_ew(j,2) = 10.d0 !vel !east
       v_bc_ew(j,2) = 10.d0
c       dp_bc_ew(j,2) = 0.d0
 
      enddo

      do i = 1, nx

        u_bc_ns(i,1) = vel !0.d0 ! south 
        v_bc_ns(i,1) = vel !0.d0 
c        dp_bc_ns(i,1) = 0.d0

        u_bc_ns(i,2) =  10.d0 ! north 
        v_bc_ns(i,2) = 10.d0 
c        dp_bc_ns(i,2) = 0.d0

      enddo

      return
      end
