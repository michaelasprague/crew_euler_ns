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

      subroutine dp_edges(dpdx_edge, dpdy_edge, p,
     &                    dp_bc_ew, dp_bc_ns, 
     &                    nx, ny, dx, dy)
      
C  Calculate dp/dx and dp/dy at cell edges
 
      implicit double precision (a-h,o-z)

C Input variables

      double precision p(nx,ny)
      double precision dp_bc_ew(ny,2), dp_bc_ns(nx,2)

C Output variables

      double precision dpdx_edge(nx+1, ny)
      double precision dpdy_edge(nx,   ny+1) 


      do j = 1,ny

        dpdx_edge(1,j) =  dp_bc_ew(j,1)
       
        do i = 2,nx
          dpdx_edge(i,j) = ( p(i,j) - p(i-1,j) ) /  dx
        enddo

        dpdx_edge(nx+1,j) = dp_bc_ew(j,2)

      enddo 

      do i = 1,nx
      
        dpdy_edge(i,1) = dp_bc_ns(i,1)
 
        do j = 2,ny

          dpdy_edge(i,j) = ( p(i,j) - p(i,j-1) ) /  dy

        enddo

        dpdy_edge(i,ny+1) = dp_bc_ns(i,2)

      enddo

      return
      end

     
