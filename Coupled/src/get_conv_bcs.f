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

      subroutine get_conv_bcs(u_bc_ew, u_bc_ns,
     &                        v_bc_ew, v_bc_ns,
     &                        u_bc_e_conv, u_bc_n_conv, 
     &                        v_bc_e_conv, v_bc_n_conv, 
     &                        u,v,
     &                        dt, nx, ny, dx, dy, k)

      implicit double precision (a-h,o-z)

      double precision u_bc_ew(ny,2), u_bc_ns(nx,2)
      double precision v_bc_ew(ny,2), v_bc_ns(nx,2)
      double precision dp_bc_ew(ny,2), dp_bc_ns(nx,2)

      double precision u_bc_e_conv(ny),  u_bc_n_conv(nx)
      double precision v_bc_e_conv(ny),  v_bc_n_conv(nx)

      double precision u(nx,ny)
      double precision v(nx,ny)

      ! bc_ew( number of elements, face 1-west,2-east, velocity component)
      ! bc_ns (number of elements, face 1-south, 2 north)

    
      if (k .eq. 1) then

         do j = 1,ny
           u_bc_e_conv(j) = u_bc_ew(j,2)
           v_bc_e_conv(j) = v_bc_ew(j,2)
         enddo

         do i = 1,nx
           u_bc_n_conv(i) = u_bc_ns(i,2)
           v_bc_n_conv(i) = v_bc_ns(i,2)
         enddo
       
      else

         do j = 1,ny !east 
           u_bc_e_conv(j) = ( dx*u_bc_e_conv(j) + 
     &                     2.d0*dt*u_bc_ew(j,2)*u(nx,j) ) 
     &           / (dx + 2.d0*dt*u_bc_ew(j,2))

           v_bc_e_conv(j) = ( dx*v_bc_e_conv(j) +
     &                     2.d0*dt*u_bc_ew(j,2)*v(nx,j) )
     &           / (dx + 2.d0*dt*u_bc_ew(j,2))
         enddo

         do i = 1,nx ! north
           u_bc_n_conv(i) = ( dy*u_bc_n_conv(i) +
     &                     2.d0*dt*u_bc_ns(i,2)*u(i,ny) )
     &           / (dy + 2.d0*dt*u_bc_ns(i,2))

           v_bc_n_conv(i) = ( dy*v_bc_n_conv(i) +
     &                     2.d0*dt*u_bc_ns(i,2)*v(i,ny) )
     &           / (dy + 2.d0*dt*u_bc_ns(i,2))
         enddo

      endif
  
      return
      end
