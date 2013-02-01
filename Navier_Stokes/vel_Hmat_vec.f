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

      subroutine vel_Hmat_vec(w, u,
     &                        ew_rc, ns_rc, 
     &                        nx, ny, dx, dy, dt, viscosity, neum,conv)

c velocity times matrix for momentum predictor eq. #4
c see notes pg. 10


      implicit double precision (a-h, o-z)

      double precision w(nx,ny), w1(nx,ny) ! inside arrays
      
      double precision u(nx,ny)
     
      double precision ew_rc(nx+1,ny), ns_rc(nx,ny+1)

      external Hbar_vec, A_vec
   
      logical neum
      logical conv


      a = (dx*dy)/dt

      call Hbar_vec(w, u, ew_rc, ns_rc,
     &              nx, ny, dx, dy, dt, viscosity)

      call A_vec(w1, u, ew_rc, ns_rc,
     &              nx, ny, dx, dy, dt, viscosity, neum,conv)


Ccschedule(static)
C!$omp parallel private(i,j) shared(w,a,u,nx,ny)
C!$omp do
      do j = 1,ny
        do i = 1,nx

          w(i,j) = a*u(i,j) - 0.5*( w(i,j) + w1(i,j) )

        enddo
      enddo
C!$omp end parallel

      return
      end

 
