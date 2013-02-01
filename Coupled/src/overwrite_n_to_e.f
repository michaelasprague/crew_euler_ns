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

      subroutine overwrite_n_to_e(u_e, v_e, r_e,
     &                            u_n, v_n, density, 
     &                            nx_n, ny_n, nx_e, ny_e,
     &                            nx_start, ny_start)


c overwrites Euler solution with NS solution 
c in case of two way partial bd coupling
c Notes: overwrite_n_to_e in Bits and Pieces

      implicit double precision (a-h,o-z)

      double precision u_e (nx_e, ny_e)
      double precision v_e (nx_e, ny_e)
      double precision r_e (nx_e, ny_e)

      double precision u_n(nx_n,ny_n)
      double precision v_n(nx_n,ny_n)


      do j = 1, ny_n 
        do i = 1, nx_n-1
     
          u_e(nx_start+i,ny_start+j-1) = 0.5d0* 
     &    ( u_n(i,j) + u_n(i+1,j) ) *
     &    sqrt( 2.d0*density / 
     & (r_e(nx_start+i-1,ny_start+j-1) + r_e(nx_start+i,ny_start+j-1)))

        enddo
      enddo

      do i = 1, nx_n 
        do j = 1, ny_n-1

          v_e(nx_start+i-1,ny_start+j) = 0.5d0* 
     &    ( v_n(i,j) + v_n(i,j+1) ) *
     &    sqrt( 2.d0*density / 
     & (r_e(nx_start+i-1,ny_start+j-1) + r_e(nx_start+i-1,ny_start+j)))

        enddo
      enddo



      return
      end
