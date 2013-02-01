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

      subroutine error_ns(u, v, r, u0, v0, r0,err_u,err_v,
     &                    rmsu,rmsu_rel,err_max,
     &                     nx, ny, mom)
 

      implicit double precision (a-h,o-z)

      double precision u  (nx, ny)
      double precision v  (nx, ny)
      double precision r  (nx, ny)
      double precision u0(nx,ny), v0(nx,ny), r0(nx,ny)
      double precision err_u(nx,ny)
      double precision err_v(nx,ny)

      logical mom

c if (mom) then mass flux is used for the error
     

c---- Euclidean norm of the error-----------

      rmsu = 0.d0 ! RMS of veloc error magnitude
      usum = 0.d0
      err_max = 0.d0
      u0_max = 0.d0




      do j = 1,ny
        do i =1,nx

         if (mom) then
          err_u(i,j) =   u(i,j)*r(i,j) - u0(i,j)*r0(i,j)
          err_v(i,j) =   v(i,j)*r(i,j) - v0(i,j)*r0(i,j) !mass flux 
          u0_mag_sq = (u0(i,j)*r0(i,j))**2 + (v0(i,j)*r0(i,j))**2    
         else
          err_u(i,j) =   u(i,j) - u0(i,j)
          err_v(i,j) =   v(i,j) - v0(i,j) ! velocity
          u0_mag_sq = u0(i,j)**2 + v0(i,j)**2    
         endif

          err_mag_sq = err_u(i,j)**2 + err_v(i,j)**2
 
          rmsu = rmsu + err_mag_sq 
          
          usum = usum + u0_mag_sq
    
          err_mag = sqrt(err_mag_sq)
          u0_mag = sqrt(u0_mag_sq)
                   
          if (err_mag .gt. err_max) err_max = err_mag
          if (u0_mag .gt. u0_max) u0_max = u0_mag
 

        enddo
      enddo

      rmsu = sqrt(rmsu / float(nx*ny))

      rmsu_rel = sqrt(rmsu/usum)
      err_max = err_max/u0_max

 
      return
      end
