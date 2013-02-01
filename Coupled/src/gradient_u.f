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

      subroutine gradient_u(du_dx,du_dy,u,dx,dy,nx,ny)

      implicit double precision (a-h,o-z)

      double precision u (nx, ny)
 
      double precision du_dx(nx-1, ny)
      double precision du_dy(nx, ny-1)


      do j = 1, ny
        do i = 1, nx-1
           
           du_dx(i,j) = ( u(i+1,j) - u(i,j) )/dx           

        enddo
      enddo

      do i = 1,nx
        do j = 1, ny-1

          du_dy(i,j) = ( u(i,j+1) - u(i,j) )/ dy

        enddo
      enddo

    

      return
      end
