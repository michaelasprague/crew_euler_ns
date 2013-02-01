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

      subroutine pres_rhs(p_rhs, ew, ns, u, v,
     &                    u_bc_ew, u_bc_ns,
     &                    v_bc_ew, v_bc_ns,
     &                    u_bc_e_conv, u_bc_n_conv,
     &                    v_bc_e_conv, v_bc_n_conv,
     &                    nx, ny, dx, dy, dt, density, neum,conv)


c calculates RHS of pressure eq.

      implicit double precision (a-h, o-z)

      double precision p_rhs(nx,ny)
      double precision u(nx,ny)
      double precision v(nx,ny)
      double precision ew(nx-1,ny), ns(nx,ny-1)
      double precision u_bc_ew(ny,2), u_bc_ns(nx,2)
      double precision v_bc_ew(ny,2), v_bc_ns(nx,2)
      double precision u_bc_e_conv(ny), u_bc_n_conv(nx)
      double precision v_bc_e_conv(ny), v_bc_n_conv(nx)

c inside arrays
      double precision u_bc_ew_rhs(ny,2), u_bc_ns_rhs(nx,2)
      double precision v_bc_ew_rhs(ny,2), v_bc_ns_rhs(nx,2)
   
      logical neum
      logical conv

      

      a   = (dx * dy) / dt
      ro  = density
      dy2 = dy * dy
      dx2 = dx * dx

c--------regular----------------
         do j = 1,ny
           u_bc_ew_rhs(j,1) = u_bc_ew(j,1)
           u_bc_ew_rhs(j,2) = u_bc_ew(j,2)
           v_bc_ew_rhs(j,1) = v_bc_ew(j,1)
           v_bc_ew_rhs(j,2) = v_bc_ew(j,2)
         enddo

         do i = 1,nx
           u_bc_ns_rhs(i,1) = u_bc_ns(i,1)
           u_bc_ns_rhs(i,2) = u_bc_ns(i,2)
           v_bc_ns_rhs(i,1) = v_bc_ns(i,1)
           v_bc_ns_rhs(i,2) = v_bc_ns(i,2)

         enddo

c---------neumann--------------
c instead of us used u_no_pres before - looked the same
c still think both are not right.
c Mike says leave it alone!

 
      if (neum) then

         do j = 1,ny
c           u_bc_ew_rhs(j,2) = 1.5d0*u(nx,j)-0.5d0*u(nx-1,j)
c           v_bc_ew_rhs(j,2) = 1.5d0*v(nx,j)-0.5d0*v(nx-1,j)
           u_bc_ew_rhs(j,2) = u(nx,j)
           v_bc_ew_rhs(j,2) = v(nx,j)
 
         enddo

         do i = 1,nx
c           u_bc_ns_rhs(i,2) = 1.5d0*u(i,ny)-0.5d0*u(i,ny-1)
c           v_bc_ns_rhs(i,2) = 1.5d0*v(i,ny)-0.5d0*v(i,ny-1) 
           u_bc_ns_rhs(i,2) = u(i,ny)
           v_bc_ns_rhs(i,2) = v(i,ny)

         enddo

c-------convective-------------- !trash!
c       else if (conv) then
c
c         do j = 1,ny !east 
c           u_bc_ew_rhs(j,2) = ( dx*u_bc_e_conv(j) +
c     &                     2*dt*u_bc_ew(j,2)*us(nx,j) )
c     &           / (dx + 2.d0*dt*u_bc_ew(j,2))
c
c           v_bc_ew_rhs(j,2) = ( dx*v_bc_e_conv(j) +
c     &                     2*dt*u_bc_ew(j,2)*vs(nx,j) )
c     &           / (dx + 2.d0*dt*u_bc_ew(j,2))
c         enddo
c
c         do i = 1,nx ! north
c           u_bc_ns_rhs(i,2) = ( dy*u_bc_n_conv(i) +
c     &                     2*dt*u_bc_ns(i,2)*us(i,ny) )
c     &           / (dy + 2.d0*dt*u_bc_ns(i,2))
c
c           v_bc_ns_rhs(i,2) = ( dy*v_bc_n_conv(i) +
c     &                     2*dt*u_bc_ns(i,2)*vs(i,ny) )
c     &           / (dy + 2.d0*dt*u_bc_ns(i,2))
c         enddo

       endif
 
 
c Inside

      do j = 2,ny-1

        do i = 2,nx-1

          p_rhs(i,j) = 
     &      + ew(i,j)   
     &      - ew(i-1,j) 
     &      + ns(i,j)   
     &      - ns(i,j-1) 

        enddo
  
      enddo


c South
      j = 1
      do i = 2, nx-1

        p_rhs(i,j) = 
     &      + ew(i,j)   
     &      - ew(i-1,j) 
     &      + ns(i,j)   
     &      - v_bc_ns_rhs(i,1) * dx 

      enddo

c North
      j = ny
      do i = 2,nx-1


        p_rhs(i,j) =
     &      + ew(i,j)   
     &      - ew(i-1,j) 
     &      + v_bc_ns_rhs(i,2) * dx 
     &      - ns(i,j-1) 
        

      enddo
  
c West

      i = 1
      do j = 2,ny-1
        p_rhs(i,j) = 
     &    + ew(i,j)   
     &    - u_bc_ew_rhs(j,1) * dy
     &    + ns(i,j)   
     &    - ns(i,j-1) 
      enddo

c East
      
      i = nx
      do j = 2,ny-1
    
        p_rhs(i,j) = 
     &    + u_bc_ew_rhs(j,2) * dy 
     &    - ew(i-1,j) 
     &    + ns(i,j)   
     &    - ns(i,j-1) 
   

      enddo
 
c Corners


c SW
      i = 1
      j = 1
       p_rhs(i,j) = 
     &            + ew(i,j)
     &            - u_bc_ew_rhs(j,1) * dy 
     &            + ns(i,j)
     &            - v_bc_ns_rhs(i,1) * dx 

c SE
      j = 1
      i = nx
 
      p_rhs(i,j) = 
     &     + u_bc_ew_rhs(j,2) * dy 
     &     - ew(i-1,j)
     &     + ns(i,j)  
     &     - v_bc_ns_rhs(i,1) * dx 



c NE    

      j = ny
      i = nx

      p_rhs(i,j) = 
     &    + u_bc_ew_rhs(j,2) * dy 
     &    - ew(i-1,j) 
     &    + v_bc_ns_rhs(i,2) * dx 
     &    - ns(i,j-1) 
     


c NW 

      j = ny
      i = 1
   
      p_rhs(i,j) = 
     &    + ew(i,j)  
     &    - u_bc_ew_rhs(j,1) * dy 
     &    + v_bc_ns_rhs(i,2) * dx 
     &    - ns(i,j-1) 



c multiply by -1 for pos. def. p_mat_vec is also multiplied by -1

      do i = 1,nx
        do j = 1,ny
          p_rhs(i,j) = -1.d0*p_rhs(i,j) 
        enddo
      enddo

c Done!

      return
      end

 
