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

      subroutine vel_rhs(rhs, u, Hnun, Anun, bn, bn1,
     &                   f, dp,
     &                   nx, ny, dx, dy, dt, k, 
     &                   density)

c calculates the rhs of momentum eq.
c use k = 1 for u; and k=2 for v;
c see pg. #7 for notes

      implicit double precision (a-h, o-z)

      double precision rhs(nx,ny), u(nx,ny)

      double precision bn (2*nx + 2*ny - 4, 2)
      double precision bn1(2*nx + 2*ny - 4, 2)

      double precision f(nx,ny,2), dp(nx,ny)


      double precision Hnun(nx,ny), Anun(nx,ny) ! inside array

      dxdy = dx * dy
      a    = dxdy / dt
 
      do j = 1,ny

        do i = 1,nx

         rhs(i,j) = 0.5d0 * ( Hnun(i,j)   + Anun(i,j) ) 
     &            + dxdy * ( f(i,j,k) - dp(i,j) / density )
     &            + a * u(i,j)

        enddo

      enddo


c now add BC's;  right-hand around the boundary

c South
      do i = 1,nx
    
        rhs(i,1) = rhs(i,1) + 0.5d0 * ( bn(i,k) + bn1(i,k) )

      enddo


c East
      do i = nx+1, nx+ny-1
   
         j = i - nx + 1 ! so j = 2,ny

         rhs(nx,j) = rhs(nx,j) + 0.5d0 * ( bn(i,k) + bn1(i,k) )
 
      enddo

c North 
      do i = nx+ny, 2*nx+ny-2

         j = -(i-2*nx-ny+1) ! so j = nx-1,1

         rhs(j,ny) = rhs(j,ny) + 0.5d0 * ( bn(i,k) + bn1(i,k) )  

      enddo

c West
      do i = 2*nx+ny-1, 2*nx+2*ny-4
   
        j = -(i-2*nx-2*ny+2) ! so j = ny-1, 2

        rhs(1,j) = rhs(1,j) + 0.5d0 * ( bn(i,k) + bn1(i,k) )

      enddo

      return
      end

