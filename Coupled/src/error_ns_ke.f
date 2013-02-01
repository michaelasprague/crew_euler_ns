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

      subroutine error_ns_ke(e, e0, err_e,
     &                       rmse,rmse_rel,
     &                       nx, ny)
 

      implicit double precision (a-h,o-z)

      double precision e  (nx, ny)
      double precision e0  (nx, ny)
      double precision err_e(nx,ny)




      rmse = 0.d0 ! RMS of veloc error magnitude
      esum = 0.d0
      err_max = 0.d0
      e_max = 0.d0




      do j = 1,ny
        do i =1,nx

          err_e(i,j) =   e(i,j) - e0(i,j)

          err_e_sq = err_e(i,j)**2 
 
          rmse = rmse + err_e_sq 
          
          esum = esum + e0(i,j)**2
    
        enddo
      enddo

      rmse = sqrt(rmse / float(nx*ny))

      rmse_rel = sqrt(rmse/esum)
 
      return
      end
