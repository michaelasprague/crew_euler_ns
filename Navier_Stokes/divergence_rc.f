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

      subroutine divergence_rc(div, u, v,
     &                      u_bc_ew, v_bc_ns,
     &                      dpdx, dpdy,
     &                      dpdx_edge, dpdy_edge, 
     &                      An, A_ew, A_ns,
     &                      div_rms, div_max, density,
     &                      dx, dy, dt, nx, ny)

      implicit double precision(a-h,o-z)

      double precision div(nx,ny), u(nx,ny), v(nx,ny)
      double precision u_bc_ew(ny,2)
      double precision v_bc_ns(nx,2)

      double precision An(nx,ny) ! A^n times ones 
      double precision A_ew(nx-1, ny  )
      double precision A_ns(nx,   ny-1)

      double precision dpdx(nx,ny), dpdy(nx,ny)
      double precision dpdx_edge(nx+1,ny), dpdy_edge(nx,ny+1)

      double precision u_ew(nx+1, ny)
      double precision v_ns(nx, ny+1)

      a = (dx*dy)/dt
   
      do j = 1, ny
        
        u_ew(1,j) = u_bc_ew(j,1)

        do i = 2, nx
           u_ew(i,j) = 
     &      0.5*(  
     &      + u(i-1,j) 
     &         + dx*dy*dpdx(i-1,j) / (density * (a- 0.5d0*An(i-1,j) ))
     &      + u(i,j)  
     &         + dx*dy*dpdx(i,j)   / (density * (a- 0.5d0*An(i  ,j) ))
     &          )
     &        - dx*dy*dpdx_edge(i,j) * A_ew(i-1,j) / density

        enddo

        u_ew(nx+1,j) = u_bc_ew(j,2)

      enddo

      do i = 1, nx
        
        v_ns(i,1) = v_bc_ns(i,1)

        do j = 2, ny
           v_ns(i,j) = 
     &      0.5*(  
     &      + v(i,j-1) 
     &         + dx*dy*dpdy(i,j-1) / (density * (a- 0.5d0*An(i,j-1) ))
     &      + v(i,j)  
     &         + dx*dy*dpdy(i,j)   / (density * (a- 0.5d0*An(i  ,j) ))
     &          )
     &        - dx*dy*dpdy_edge(i,j) * A_ns(i,j-1) / density 

        enddo

        v_ns(i,ny+1) = v_bc_ns(i,2)

      enddo

      do j = 1, ny
        do i = 1, nx
          div(i,j) = dy*(u_ew(i,j) - u_ew(i+1,j))
     &             + dx*(v_ns(i,j) - v_ns(i,j+1))
        enddo
      enddo

      div_max = 0.d0
      div_rms = 0.d0
      do j = 1, ny
        do i = 1, nx
          div(i,j) = abs(div(i,j))
          if (div(i,j) .gt. div_max) div_max = div(i,j)
          div_rms = div_rms + div(i,j)**2
        enddo
      enddo

      div_rms = sqrt(div_rms / float(nx*ny))

      return
      end
