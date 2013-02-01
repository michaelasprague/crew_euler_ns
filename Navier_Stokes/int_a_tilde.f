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

      subroutine int_a_tilde(A_ew, A_ns, AH, dx, dy, dt, nx, ny,
     &                       A_e, A_n)

c interpolates to the inside edges (a_tilda)^-1 
c needed for p_mat_vec

      implicit double precision (a-h,o-z)
 
      double precision AH(nx, ny)
      double precision w(nx, ny)
      double precision A_ew(nx-1, ny  )
      double precision A_ns(nx,   ny-1)
      double precision A_e(ny) ! (a_tilda_p)^-1 along east edge 
      double precision A_n(nx)! ! (a_tilda_p)^-1 along north edge

      a = (dx*dy)/dt

c calculating (a_tilde)^-1

      do j = 1,ny
        do i = 1,nx

          w(i,j) = 1.d0 / (a - 0.5d0*AH(i,j) )

        enddo
      enddo


c interpolating to e-w edges

      do j = 1, ny
        do i = 1, nx-1

          A_ew(i,j) = 0.5d0*( w(i,j) + w(i+1,j) )

        enddo
      enddo

c interpolating to the n-s edges


      do j = 1, ny-1
        do i = 1, nx

          A_ns(i,j) = 0.5d0*( w(i,j) + w(i,j+1) )

        enddo
      enddo

c taking w only along East and North edge

      do i = 1, nx
        A_n(i) = w(i,ny)
      enddo

      do j = 1, ny
        A_e(j) = w(nx,j)
      enddo


   
      return
      end
         
