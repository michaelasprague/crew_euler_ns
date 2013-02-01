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

      subroutine rhie_chow(u, v, p, ew_rc, ns_rc,
     &                 dpdx, dpdy,
     &                 dpdx_edge, dpdy_edge, 
     &                 An, A_ew, A_ns,
     &                 u_bc_ew, u_bc_ns, v_bc_ew, v_bc_ns,
     &                 dx, dy, dt, density, nx, ny, rc, neum)

      implicit double precision (a-h,o-z)

C Rhie-Chow interpolated velocites PLUS known boundary values

C primary variables
      double precision u(nx,ny), v(nx,ny), p(nx,ny)

      double precision dpdx(nx,ny)
      double precision dpdy(nx,ny)

      double precision dpdx_edge(nx+1,ny)
      double precision dpdy_edge(nx,ny+1)


C   Storage for Rhi-Chow interpolated velocities
      double precision ew_rc  (nx+1, ny)
      double precision ns_rc  (nx  , ny+1)

      double precision An(nx,ny) ! A^* times ones
      double precision A_ew(nx-1, ny  )
      double precision A_ns(nx,   ny-1)

C   boundary condition stuff
      double precision u_bc_ew(ny,2), v_bc_ew(ny,2)
      double precision u_bc_ns(nx,2), v_bc_ns(nx,2)

      logical rc
      logical neum

      a = (dx*dy)/dt

C Vertical sides including boundaries
      do j = 1,ny

        do i = 2,nx

          if (rc) then 
            ew_rc(i,j) = 
     &        0.5*(
     &        + u(i-1,j)
     &           + dx*dy*dpdx(i-1,j) / (density * (a- 0.5d0*An(i-1,j) ))
     &        + u(i,j)
     &           + dx*dy*dpdx(i,j)   / (density * (a- 0.5d0*An(i  ,j) ))
     &            )
     &          - dx*dy*dpdx_edge(i,j) * A_ew(i-1,j) / density
          else
            ew_rc(i,j) = 
     &        0.5*( + u(i-1,j) + u(i,j) )
          endif


        enddo

        ew_rc(1,   j) = u_bc_ew(j,1) !west

        if (neum) then
          ew_rc(nx+1,j) = u(nx,j) !east
        else 
          ew_rc(nx+1,j) = u_bc_ew(j,2) !east
        endif
 
      enddo

C horizontal sides including boundaries
      do i = 1,nx

        do j = 2,ny


          if (rc) then 
            ns_rc(i,j) = 
     &        0.5*(
     &        + v(i,j-1)
     &           + dx*dy*dpdy(i,j-1) / (density * (a- 0.5d0*An(i,j-1) ))
     &        + v(i,j)
     &           + dx*dy*dpdy(i,j)   / (density * (a- 0.5d0*An(i  ,j) ))
     &            )
     &          - dx*dy*dpdy_edge(i,j) * A_ns(i,j-1) / density
           else
            ns_rc(i,j) = 
     &        0.5*( + v(i,j-1) + v(i,j) )
       
           endif
        enddo


        ns_rc(i,   1) = v_bc_ns(i,1)  ! south
        
        if (neum) then
          ns_rc(i,ny+1) = v(i,ny)  ! north
        else 
          ns_rc(i,ny+1) = v_bc_ns(i,2)  ! north
        endif

      enddo

      return
      end
