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

      subroutine add_ke( u_bc_ew, u_bc_ns, v_bc_ew, v_bc_ns,
     &                   ue_bc_ew, ue_bc_ns, ve_bc_ew, ve_bc_ns,
     &                   re_bc_ns, re_bc_ew,
     &                   nx_n,ny_n,dx_n,dy_n,density)

c calculates constat vector (u,v) to add to NS boundaries
c so that KE is the same as in Euler


      implicit double precision (a-h,o-z)


      double precision u_bc_ew(ny_n,2), u_bc_ns(nx_n,2)
      double precision v_bc_ew(ny_n,2), v_bc_ns(nx_n,2)

      double precision re_bc_ew(ny_n,2), re_bc_ns(nx_n,2)
      double precision ue_bc_ew(ny_n,2), ue_bc_ns(nx_n,2)
      double precision ve_bc_ew(ny_n,2), ve_bc_ns(nx_n,2)


      closed_integral_u = 0.d0
      closed_integral_v = 0.d0

c length of the NS boundary
      gamma_length = 2.d0*nx_n*dx_n + 2.d0*ny_n*dy_n

       

      do i = 1,nx_n

        ! u velocity 
        closed_integral_u = closed_integral_u +
     &    ( sqrt(re_bc_ns(i,1))*ue_bc_ns(i,1) -   !south
     &      sqrt(density)*u_bc_ns(i,1) ) * dx_n +
     &    ( sqrt(re_bc_ns(i,2))*ue_bc_ns(i,2) -   !north
     &      sqrt(density)*u_bc_ns(i,2) ) * dx_n

        ! v velocity
        closed_integral_v = closed_integral_v +
     &    ( sqrt(re_bc_ns(i,1))*ve_bc_ns(i,1) -   !south
     &      sqrt(density)*v_bc_ns(i,1) ) * dx_n +
     &    ( sqrt(re_bc_ns(i,2))*ve_bc_ns(i,2) -   !north
     &      sqrt(density)*v_bc_ns(i,2) ) * dx_n

 
      enddo



      do j = 1,ny_n

        ! u velocity 
        closed_integral_u = closed_integral_u +
     &    ( sqrt(re_bc_ew(j,1))*ue_bc_ew(j,1) -   !west
     &      sqrt(density)*u_bc_ew(j,1) ) * dy_n +
     &    ( sqrt(re_bc_ew(j,2))*ue_bc_ew(j,2) -   !east
     &      sqrt(density)*u_bc_ew(j,2) ) * dy_n

        ! v velocity
        closed_integral_v = closed_integral_v +
     &    ( sqrt(re_bc_ew(j,1))*ve_bc_ew(j,1) -   !east
     &      sqrt(density)*v_bc_ew(j,1) ) * dy_n +
     &    ( sqrt(re_bc_ew(j,2))*ve_bc_ew(j,2) -   !west
     &      sqrt(density)*v_bc_ew(j,2) ) * dy_n


      enddo


c calc constant u and v
      u = closed_integral_u/gamma_length
      v = closed_integral_v/gamma_length



c add u and v to NS boundaries


      do i = 1,nx_n

        u_bc_ns(i,1) = u_bc_ns(i,1) + u !south
        u_bc_ns(i,2) = u_bc_ns(i,2) + u !north

        v_bc_ns(i,1) = v_bc_ns(i,1) + v !south
        v_bc_ns(i,2) = v_bc_ns(i,2) + v !north

      enddo

      do j = 1, ny_n

        u_bc_ew(j,1) = u_bc_ew(j,1) + u !west
        u_bc_ew(j,2) = u_bc_ew(j,2) + u !east
      
        v_bc_ew(j,1) = v_bc_ew(j,1) + v !west
        v_bc_ew(j,2) = v_bc_ew(j,2) + v !east

      enddo

        
      return
      end



