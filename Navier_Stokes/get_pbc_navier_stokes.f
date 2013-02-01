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

      subroutine get_pbc_navier_stokes(u,v,
     &                   u_bc_ew_n, u_bc_ns_n, v_bc_ew_n, v_bc_ns_n,
     &                   u_bc_ew_n1, u_bc_ns_n1, v_bc_ew_n1, v_bc_ns_n1,
     &                   dp_bc_ew, dp_bc_ns, 
     &                   dt, nx, ny, 
     &                   dx, dy, density, viscosity, zero_dpdn)

! calculates pressure-gradient bc that satisfies the navier-stokes equations

      implicit double precision (a-h,o-z)

      double precision u(nx, ny)
      double precision v(nx, ny)

      double precision u_bc_ew_n(ny,2), u_bc_ns_n(nx,2)
      double precision v_bc_ew_n(ny,2), v_bc_ns_n(nx,2)

      double precision u_bc_ew_n1(ny,2), u_bc_ns_n1(nx,2)
      double precision v_bc_ew_n1(ny,2), v_bc_ns_n1(nx,2)

      double precision dp_bc_ew(ny,2), dp_bc_ns(nx,2)

      logical zero_dpdn


      if (zero_dpdn) then
       do j = 1, 2
         do i = 1, nx
           dp_bc_ns(i,j) = 0.d0
         enddo
         do i = 1, ny
           dp_bc_ew(i,j) = 0.d0
         enddo
       enddo

      else 


! West & East sides
      do i = 2, ny - 1


! --------------- West -----------------
        dudt = (u_bc_ew_n1(i,1) - u_bc_ew_n(i,1)) / dt

        dudx = (u(1,i) - u_bc_ew_n1(i,1) ) / (0.5*dx)

        dudy = (u_bc_ew_n1(i+1,1) - u_bc_ew_n1(i-1,1)  ) / (2.d0 * dy)

        ududx = u_bc_ew_n1(i,1) * dudx

        vdudy = v_bc_ew_n1(i,1) * dudy

        d2udx2 = (+ 2.d0 * u_bc_ew_n1(i,1)
     &            - 3.d0 * u(1,i)
     &             + 1.d0 * u(2,i)) * 4.d0 / (3.d0 * dx * dx)

        d2udy2 = (u_bc_ew_n1(i+1,1) 
     &            - 2.d0 * u_bc_ew_n1(i,1)
     &            +  u_bc_ew_n1(i-1,1) ) / (dy*dy)
 
        dp_bc_ew(i,1) = density*(-dudt - ududx - vdudy 
     &                           + viscosity * (d2udx2 + d2udy2) )

 
! --------------- East -----------------
        dudt =  (u_bc_ew_n1(i,2) - u_bc_ew_n(i,2)) / dt

        dudx = -(u(nx,i) - u_bc_ew_n1(i,2) ) / (0.5*dx)

        dudy = (u_bc_ew_n1(i+1,2) - u_bc_ew_n1(i-1,2)  ) / (2.d0 * dy)

        ududx = u_bc_ew_n1(i,2) * dudx

        vdudy = v_bc_ew_n1(i,2) * dudy

        d2udx2 = (+ 2.d0 * u_bc_ew_n1(i,2)
     &            - 3.d0 * u(nx,i)
     &            + 1.d0 * u(nx-1,i)) * 4.d0 / (3.d0 * dx * dx)

        d2udy2 = (u_bc_ew_n1(i+1,2) 
     &            - 2.d0 * u_bc_ew_n1(i,2)
     &            +  u_bc_ew_n1(i-1,2) ) / (dy*dy)
 
        tmp = dp_bc_ew(i,2) 
        dp_bc_ew(i,2) = density*(-dudt - ududx - vdudy 
     &                           + viscosity * (d2udx2 + d2udy2) )



      enddo


!----------------------------------------------------------------------------

! North & South Sides
      do i = 2, nx-1

! --------------- South -----------------

        dvdt =  (v_bc_ns_n1(i,1) - v_bc_ns_n(i,1)) / dt

        dvdx = (v_bc_ns_n1(i+1,1) - v_bc_ns_n1(i-1,1)  ) / (2.d0 * dx)

        dvdy = (v(i,1) - v_bc_ns_n1(i,1) ) / (0.5*dy)

        udvdx = u_bc_ns_n1(i,1) * dvdx

        vdvdy = v_bc_ns_n1(i,1) * dvdy

        d2vdy2 = (+ 2.d0 * v_bc_ns_n1(i,1)
     &            - 3.d0 * v(i,1)
     &            + 1.d0 * v(i,2)) * 4.d0 / (3.d0 * dy * dy)

        d2vdx2 = (v_bc_ns_n1(i+1,1) 
     &            - 2.d0 * v_bc_ns_n1(i,1)
     &            + v_bc_ns_n1(i-1,1) ) / (dx*dx)
 
        dp_bc_ns(i,1) = density*(-dvdt - udvdx - vdvdy 
     &                           + viscosity * (d2vdx2 + d2vdy2) )
 
! --------------- North -----------------
        dvdt =  (v_bc_ns_n1(i,2) - v_bc_ns_n(i,2)) / dt

        dvdx = (v_bc_ns_n1(i+1,2) - v_bc_ns_n1(i-1,2)  ) / (2.d0 * dx)

        dvdy = (v(i,ny) - v_bc_ns_n1(i,2) ) / (0.5*dy)

        udvdx = u_bc_ns_n1(i,2) * dvdx

        vdvdy = v_bc_ns_n1(i,2) * dvdy

        d2vdy2 = (+ 2.d0 * v_bc_ns_n1(i,2)
     &            - 3.d0 * v(i,ny)
     &            + 1.d0 * v(i,ny-1)) * 4.d0 / (3.d0 * dy * dy)

        d2vdx2 = (v_bc_ns_n1(i+1,2) 
     &            - 2.d0 * v_bc_ns_n1(i,2)
     &            + v_bc_ns_n1(i-1,2) ) / (dx*dx)
 
        tmp = dp_bc_ns(i,2) 
        dp_bc_ns(i,2) = density*(-dvdt - udvdx - vdvdy 
     &                           + viscosity * (d2vdx2 + d2vdy2) )


      enddo

! ---------------------------------------------------------------------
! West-South Corner
        i = 1
        dudt = (u_bc_ew_n1(1,1) - u_bc_ew_n(1,1)) / dt

        dudx = (u(1,1) - u_bc_ew_n1(1,1) ) / (0.5*dx)

        dudy = (u_bc_ew_n1(2,1) - u_bc_ew_n1(1,1)  ) / dy

        ududx = u_bc_ew_n1(1,1) * dudx

        vdudy = v_bc_ew_n1(1,1) * dudy

        d2udx2 = (+ 2.d0 * u_bc_ew_n1(1,1)
     &            - 3.d0 * u(1,1)
     &             + 1.d0 * u(2,1)) * 4.d0 / (3.d0 * dx * dx)

        d2udy2 = (u_bc_ew_n1(3,1) 
     &            - 2.d0 * u_bc_ew_n1(2,1)
     &            +  u_bc_ew_n1(1,1) ) / (dy*dy)
 
        dp_bc_ew(1,1) = density*(-dudt - ududx - vdudy 
     &                           + viscosity * (d2udx2 + d2udy2) )
c ------ south face
        dvdt =  (v_bc_ns_n1(1,1) - v_bc_ns_n(1,1)) / dt

        dvdx = (v_bc_ns_n1(2,1) - v_bc_ns_n1(1,1)  ) / dx

        dvdy = (v(1,1) - v_bc_ns_n1(1,1) ) / (0.5*dy)

        udvdx = u_bc_ns_n1(1,1) * dvdx

        vdvdy = v_bc_ns_n1(1,1) * dvdy

        d2vdy2 = (+ 2.d0 * v_bc_ns_n1(1,1)
     &            - 3.d0 * v(1,1)
     &            + 1.d0 * v(1,2)) * 4.d0 / (3.d0 * dy * dy)

        d2vdx2 = (v_bc_ns_n1(3,1) 
     &            - 2.d0 * v_bc_ns_n1(2,1)
     &            + v_bc_ns_n1(3,1) ) / (dx*dx)
 
        dp_bc_ns(1,1) = density*(-dvdt - udvdx - vdvdy 
     &                           + viscosity * (d2vdx2 + d2vdy2) )
        
! ---------------------------------------------------------------------
! West-North Corner
! ---- west face
        dudt = (u_bc_ew_n1(ny,1) - u_bc_ew_n(ny,1)) / dt

        dudx = (u(1,ny) - u_bc_ew_n1(ny,1) ) / (0.5*dx)

        dudy = (u_bc_ew_n1(ny,1) - u_bc_ew_n1(ny-1,1)  ) / dy

        ududx = u_bc_ew_n1(ny,1) * dudx

        vdudy = v_bc_ew_n1(ny,1) * dudy

        d2udx2 = (+ 2.d0 * u_bc_ew_n1(ny,1)
     &            - 3.d0 * u(1,ny)
     &             + 1.d0 * u(2,ny)) * 4.d0 / (3.d0 * dx * dx)

        d2udy2 = (u_bc_ew_n1(ny,1) 
     &            - 2.d0 * u_bc_ew_n1(ny-1,1)
     &            +  u_bc_ew_n1(ny-2,1) ) / (dy*dy)
 
        dp_bc_ew(1,1) = density*(-dudt - ududx - vdudy 
     &                           + viscosity * (d2udx2 + d2udy2) )

! ---- north face
        dvdt =  (v_bc_ns_n1(1,2) - v_bc_ns_n(1,2)) / dt

        dvdx = (v_bc_ns_n1(2,2) - v_bc_ns_n1(1,2)  ) / dx

        dvdy = -(v(1,ny) - v_bc_ns_n1(1,2) ) / (0.5*dy)

        udvdx = u_bc_ns_n1(1,2) * dvdx

        vdvdy = v_bc_ns_n1(1,2) * dvdy

        d2vdy2 = (+ 2.d0 * v_bc_ns_n1(1,2)
     &            - 3.d0 * v(1,ny)
     &            + 1.d0 * v(1,ny-1)) * 4.d0 / (3.d0 * dy * dy)

        d2vdx2 = (v_bc_ns_n1(3,2) 
     &            - 2.d0 * v_bc_ns_n1(2,2)
     &            + v_bc_ns_n1(1,2) ) / (dx*dx)
 
        dp_bc_ns(1,2) = density*(-dvdt - udvdx - vdvdy 
     &                           + viscosity * (d2vdx2 + d2vdy2) )

! ---------------------------------------------------------------------
! East-North Corner
! ---- east face
        dudt = (u_bc_ew_n1(ny,2) - u_bc_ew_n(ny,2)) / dt

        dudx = -(u(nx,ny) - u_bc_ew_n1(ny,2) ) / (0.5*dx)

        dudy = (u_bc_ew_n1(ny,2) - u_bc_ew_n1(ny-1,2)  ) / dy

        ududx = u_bc_ew_n1(ny,2) * dudx

        vdudy = v_bc_ew_n1(ny,2) * dudy

        d2udx2 = (+ 2.d0 * u_bc_ew_n1(ny,2)
     &            - 3.d0 * u(nx,ny)
     &             + 1.d0 * u(nx-1,ny)) * 4.d0 / (3.d0 * dx * dx)

        d2udy2 = (u_bc_ew_n1(ny,2) 
     &            - 2.d0 * u_bc_ew_n1(ny-1,2)
     &            +  u_bc_ew_n1(ny-2,2) ) / (dy*dy)
 
        dp_bc_ew(ny,2) = density*(-dudt - ududx - vdudy 
     &                           + viscosity * (d2udx2 + d2udy2) )

! ---- north face
        dvdt =  (v_bc_ns_n1(nx,2) - v_bc_ns_n(nx,2)) / dt

        dvdx = (v_bc_ns_n1(nx,2) - v_bc_ns_n1(nx-1,2)  ) / dx

        dvdy = -(v(nx,ny) - v_bc_ns_n1(nx,2) ) / (0.5*dy)

        udvdx = u_bc_ns_n1(nx,2) * dvdx

        vdvdy = v_bc_ns_n1(nx,2) * dvdy

        d2vdy2 = (+ 2.d0 * v_bc_ns_n1(nx,2)
     &            - 3.d0 * v(nx,ny)
     &            + 1.d0 * v(nx,ny-1)) * 4.d0 / (3.d0 * dy * dy)

        d2vdx2 = (v_bc_ns_n1(nx,2) 
     &            - 2.d0 * v_bc_ns_n1(nx-1,2)
     &            + v_bc_ns_n1(nx-2,2) ) / (dx*dx)
 
        dp_bc_ns(nx,2) = density*(-dvdt - udvdx - vdvdy 
     &                           + viscosity * (d2vdx2 + d2vdy2) )

! ---------------------------------------------------------------------
! East-South Corner
! ---- east face
        dudt = (u_bc_ew_n1(1,2) - u_bc_ew_n(1,2)) / dt

        dudx = -(u(nx,1) - u_bc_ew_n1(1,2) ) / (0.5*dx)

        dudy = -(u_bc_ew_n1(1,2) - u_bc_ew_n1(2,2)  ) / dy

        ududx = u_bc_ew_n1(1,2) * dudx

        vdudy = v_bc_ew_n1(1,2) * dudy

        d2udx2 = (+ 2.d0 * u_bc_ew_n1(1,2)
     &            - 3.d0 * u(nx,1)
     &             + 1.d0 * u(nx-1,1)) * 4.d0 / (3.d0 * dx * dx)

        d2udy2 = (u_bc_ew_n1(1,2) 
     &            - 2.d0 * u_bc_ew_n1(2,2)
     &            +  u_bc_ew_n1(3,2) ) / (dy*dy)
 
        dp_bc_ew(1,2) = density*(-dudt - ududx - vdudy 
     &                           + viscosity * (d2udx2 + d2udy2) )

! ---- south face
        dvdt =  (v_bc_ns_n1(nx,1) - v_bc_ns_n(nx,1)) / dt

        dvdx = (v_bc_ns_n1(nx,1) - v_bc_ns_n1(nx-1,1)  ) / dx

        dvdy =  (v(nx,1) - v_bc_ns_n1(nx,1) ) / (0.5*dy)

        udvdx = u_bc_ns_n1(nx,1) * dvdx

        vdvdy = v_bc_ns_n1(nx,1) * dvdy

        d2vdy2 = (+ 2.d0 * v_bc_ns_n1(nx,1)
     &            - 3.d0 * v(nx,1)
     &            + 1.d0 * v(nx,2)) * 4.d0 / (3.d0 * dy * dy)

        d2vdx2 = (v_bc_ns_n1(nx,1) 
     &            - 2.d0 * v_bc_ns_n1(nx-1,1)
     &            + v_bc_ns_n1(nx-2,1) ) / (dx*dx)
 
        dp_bc_ns(nx,1) = density*(-dvdt - udvdx - vdvdy 
     &                           + viscosity * (d2vdx2 + d2vdy2) )

      endif


      return
      end
