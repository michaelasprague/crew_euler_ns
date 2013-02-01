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

      subroutine init_navier_stokes(u, v, p, u_proj, v_proj,
     &         dx, dy, nx, ny, n_p)

c interpolate cell-face values from projected euler solution onto cell
c centers

      implicit double precision (a-h,o-z)

      double precision u (nx, ny)
      double precision v (nx, ny)
      double precision p (nx, ny)

      double precision u_proj (nx+2*n_p+1, ny+2*n_p)
      double precision v_proj (nx+2*n_p, ny+2*n_p+1)

c      do j = 1, ny
c        do i = 1, nx
c          u(i,j) = 0.5*(u_proj(i,j) + u_proj(i+1,j))
c          v(i,j) = 0.5*(v_proj(i,j) + v_proj(i,j+1))
c          p(i,j) = 0.d0
c        enddo 
c      enddo 

      do j = 1, ny
        do i = 1, nx
          u(i,j) = 0.5*(u_proj(i+n_p,j+n_p) + u_proj(i+n_p+1,j+n_p))
          v(i,j) = 0.5*(v_proj(i+n_p,j+n_p) + v_proj(i+n_p,j+n_p+1))
          p(i,j) = 0.d0
        enddo
      enddo


      return
      end

