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

      subroutine matlabvtk_bound(u_bc_ew, u_bc_ns,
     &                           v_bc_ew, v_bc_ns,
     &                           nx, ny, filename)

      implicit double precision (a-h,o-z)

      double precision u_bc_ew(ny,2)
      double precision u_bc_ns(nx,2)
      double precision v_bc_ew(ny,2)
      double precision v_bc_ns(nx,2)

      character filename*(*)

      open (unit=107, file=filename//'.dat', status='unknown')

      write(107,*) nx
      write(107,*) ny
     
c     u_west 
      do i = 1,nx
          write(107,*) u_bc_ew(i,1)
      enddo

c     u_east 
      do i = 1,nx
          write(107,*) u_bc_ew(i,2)
      enddo                                                         

c     u_south
      do i = 1,ny
          write(107,*) u_bc_ns(i,1)
      enddo  

c     u_north
      do i = 1,nx
          write(107,*) u_bc_ns(i,2)
      enddo

c     v_west 
      do i = 1,nx
          write(107,*) v_bc_ew(i,1)
      enddo

c     v_east 
      do i = 1,nx
          write(107,*) v_bc_ew(i,2)
      enddo                                                         

c     v_south
      do i = 1,ny
          write(107,*) v_bc_ns(i,1)
      enddo  

c     v_north
      do i = 1,nx
          write(107,*) v_bc_ns(i,2)
      enddo


c      call flush(20)

      return
      end
