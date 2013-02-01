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

      subroutine calc_error(error, error_p, 
     &            u, v, p,
     &            tu, tv, tp, 
     &            u_err_max, u_err_min,
     &            rmsu, rmsp,
     &            nx,ny)

      implicit double precision (a-h,o-z)

      double precision u  (nx, ny)
      double precision v  (nx, ny)
      double precision p  (nx, ny)
      double precision tu(nx,ny), tv(nx,ny), tp(nx,ny)
      double precision error(nx,ny), error_p(nx,ny)  ! velc and pressure error

c---- Euclidean norm of the error-----------

      rmsu = 0.d0
      rmsp = 0.d0

      u_err_max = -1d6
      u_err_min = +1d6

      do j = 1,ny
        do i =1,nx

          error(i,j) = ( ( (u(i,j) - tu(i,j)) )**2
     &                  +( (v(i,j) - tv(i,j)) )**2 ) ! velocity

          error_p(i,j) =   (p(i,j) - tp(i,j)) **2

          u_error = u(i,j) - tu(i,j) 
           
          if (u_error .gt. u_err_max)
     &          u_err_max = u_error

          if (u_error .lt. u_err_min)
     &          u_err_min = u_error

          rmsu = rmsu + error(i,j) 
          rmsp = rmsp + error_p(i,j) 

        enddo
      enddo

      if (rmsu .lt. 0) then
        write(*,*) 'rmsu less than zero; rmsu = ', rmsu
        rmsu = 0.d0
      endif
      if (rmsp .lt. 0) then
        write(*,*) 'rmsp less than zero; rmsp = ', rmsp
        rmsp = 0.d0
      endif

      rmsu = sqrt(rmsu / float(nx*ny))
      rmsp = sqrt(rmsp / float(nx*ny))

      do j = 1,ny
        do i =1,nx

          error(i,j) = sqrt(error(i,j))

        enddo
      enddo
 
      return
      end
