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

      subroutine calc_b(b, 
     &                  u_bc_ew, u_bc_ns,
     &                  v_bc_ew, v_bc_ns,
     &                  u_bc_e_conv, u_bc_n_conv,
     &                  v_bc_e_conv, v_bc_n_conv,
     &                  nx, ny, dx, dy, viscosity, neum, conv)


C calculates b terms in momentum predictor eq.#4 
c for b_n use bc's at t_n; for b_(n+1) use bc's at t_(n+1)
c See pg.8 and pg. 4-5 for notes

      implicit double precision (a-h,o-z)
  
      double precision b( 2*nx + 2*ny - 4, 2)
c  second coordinate defines velocity component: 1 for u; 2 for v

C   boundary condition stuff

      double precision u_bc_ew(ny,2), v_bc_ew(ny,2)
      double precision u_bc_ns(nx,2), v_bc_ns(nx,2)

      double precision u_bc_e_conv(ny),  u_bc_n_conv(nx)
      double precision v_bc_e_conv(ny),  v_bc_n_conv(nx)

      logical neum
      logical conv

      vi = viscosity

      do i = 1,2*nx+2*ny-4

        b(i,1) = 0.d0
        b(i,2) = 0.d0

      enddo

c South
      do i = 1, nx


        b(i,1) = b(i,1) 
     &          + u_bc_ns(i,1) * ( v_bc_ns(i,1)*dx + 2.d0 * vi*dx/dy )

        b(i,2) = b(i,2) 
     &          + v_bc_ns(i,1) * ( v_bc_ns(i,1)*dx + 2.d0 * vi*dx/dy )

      enddo


c East

      if (neum) then

        do i = nx, nx+ny-1
          b(i,1) = b(i,1) + 0.d0
          b(i,2) = b(i,2) + 0.d0
        enddo

      else if (conv) then 
        
        do i = nx, nx+ny-1

          j = i-nx+1 ! so j = 1,ny

          b(i,1) = b(i,1)
     &            + u_bc_e_conv(j) *
     &            (   (2.d0*dy*vi - dx*dy*u_bc_ew(j,2))/
     &                (dx + 2.d0*dt*u_bc_ew(j,2))        )

          b(i,2) = b(i,2)
     &             + v_bc_e_conv(j) * 
     &            (   (2.d0*dy*vi - dx*dy*u_bc_ew(j,2))/
     &                (dx + 2.d0*dt*u_bc_ew(j,2))        )

        enddo

      else
     
        do i = nx, nx+ny-1

          j = i-nx+1 ! so j = 1,ny
 
          b(i,1) = b(i,1)
     &            + u_bc_ew(j,2) * ( -u_bc_ew(j,2)*dy + 2.d0 * vi*dy/dx)

          b(i,2) = b(i,2)
     &            + v_bc_ew(j,2) * ( -u_bc_ew(j,2)*dy + 2.d0 * vi*dy/dx)
 
        enddo

      endif

c North 

      if (neum) then

      do i =   nx+ny-1, 2*nx+ny-2
        b(i,1) = b(i,1) + 0.d0
        b(i,2) = b(i,2) + 0.d0
      enddo
c missing "elseif (conv) then ...

      else

      do i =   nx+ny-1, 2*nx+ny-2

        j = -(i-2*nx-ny+1) ! so j = nx,1
        if (j .le. 0) stop 'check calc_b - spot 1'

        b(i,1) = b(i,1)
     &          + u_bc_ns(j,2) * ( -v_bc_ns(j,2)*dx + 2.d0 * vi*dx/dy)

        b(i,2) = b(i,2)
     &          + v_bc_ns(j,2) * ( -v_bc_ns(j,2)*dx + 2.d0 * vi*dx/dy)

      enddo
 
      endif



c West

      do i = 2*nx+ny-2, 2*nx+2*ny-4
 
        j = -(i-2*nx-2*ny+2) ! so j = ny,2

        if (j .le. 0) stop 'check calc_b - spot 2'

        b(i,1) = b(i,1)
     &          + u_bc_ew(j,1) * ( u_bc_ew(j,1)*dy + 2.d0*vi*dy/dx )

        b(i,2) = b(i,2)
     &          + v_bc_ew(j,1) * ( u_bc_ew(j,1)*dy + 2.d0*vi*dy/dx )
      

      enddo

c Lastly, contribution of west side to b(1,*)


      b(1,1) = b(1,1)
     &        + u_bc_ew(1,1) * ( u_bc_ew(1,1)*dy + 2.d0 * vi*dy/dx )

      b(1,2) = b(1,2) 
     &        + v_bc_ew(1,1) * ( u_bc_ew(1,1)*dy + 2.d0 * vi*dy/dx )


      return
      end
