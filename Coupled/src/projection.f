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

      subroutine projection(u,v,p,dx,dy,
     &                      nx,ny,forward, mom, density,
     &                      lag_mult_dx, lag_mult_dy, tol_cg_proj)

      implicit double precision (a-h,o-z)

      double precision u  (nx+1, ny)
      double precision v  (nx, ny+1)

      double precision p(nx, ny)  ! lagrangian multiplier
      double precision lag_mult_dx(nx+1,ny)  ! grad of lagrangian multiplier
      double precision lag_mult_dy(nx,ny+1)  ! grad of lagrangian multiplier

      double precision rhs(nx,ny)

      external proj_cg_solve, proj_mat_vec

      logical forward
      logical mom 
      logical mom_sqrt 

      cg_tol = tol_cg_proj 

      rhs_avg = 0.d0

      nx1 = nx - 1 
      ny1 = ny - 1 

      if (mom .and. .not. forward) then

        do i = 1, nx+1
           do j = 1, ny
             u(i,j) = u(i,j)*density
           enddo
         enddo

         do i = 1,nx
           do j = 1,ny+1
             v(i,j) = v(i,j)*density
           enddo
         enddo

      endif ! mom and not forward

      if (forward) then
        direction = 1.d0
      else
        direction = -1.d0  ! negative for reverse
      endif

      if (forward) then
        do j = 1, ny
          do i = 1, nx

            rhs(i,j) = (u(i+1,j)-u(i,j))/dx + (v(i,j+1) - v(i,j))/dy

            rhs_avg = rhs_avg + rhs(i,j)

            p(i,j) = 0.d0

          enddo
        enddo

        rhs_avg = rhs_avg / float ( (ny)*(nx) )

        call proj_cg_solve(p,rhs,dx,dy,nx,ny,proj_mat_vec,ierr,
     &              kfinal,cg_tol)

      endif

 
c perform projection  (forward or reverse)
      do j = 1, ny

        i = 1
! west
        dpdx = (2.*p(i,j) )/dx

        u(i,j) = u(i,j) - direction * dpdx

        lag_mult_dx(i,j) = dpdx

        do i = 2, nx
      
          dpdx = ( p(i,j) - p(i-1,j) ) / dx

          u(i,j) = u(i,j) - direction * dpdx

          lag_mult_dx(i,j) = dpdx

        enddo

        i = nx+1
! east 
        dpdx = (-2.*p(i-1,j) )/dx

        u(i,j) = u(i,j) - direction * dpdx

        lag_mult_dx(i,j) = dpdx

      enddo

      do i = 1, nx

        j = 1

! south 
!       dpdy = -(2.*p(i,j) - 3.*p(i,j+1) + p(i,j+2))/dy
!       dpdy = (p(i,j+1) - p(i,j))/dy
        dpdy = (2.*p(i,j) )/dy
c        write(*,*) 'dpdy',dpdy

        v(i,j) = v(i,j) - direction * dpdy

        lag_mult_dy(i,j) = dpdy

        do j = 2, ny
      
          dpdy = ( p(i,j) - p(i,j-1) ) / dy

          v(i,j) = v(i,j) - direction * dpdy

          lag_mult_dy(i,j) = dpdy

        enddo

        j = ny+1

! north 
!       dpdy = +(2.*p(i,j-1) - 3.*p(i,j-2) + p(i,j-3))/dy
!       dpdy = +(p(i,j) - p(i,j-1))/dy
        dpdy = (-2.*p(i,j-1) )/dy

        v(i,j) = v(i,j) - direction * dpdy
       
        lag_mult_dy(i,j) = dpdy

      enddo

      if (mom .and. forward) then
 
        do i = 1, nx+1
          do j = 1, ny
            u(i,j) = u(i,j)/density
          enddo
        enddo

        do i = 1,nx
          do j = 1,ny+1
            v(i,j) = v(i,j)/density
          enddo
        enddo

      endif ! mom and not forward


      return
      end
