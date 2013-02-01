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

      subroutine add_ke3( u_bc_ew, u_bc_ns, v_bc_ew, v_bc_ns,
     &                   ue_bc_ew, ue_bc_ns, ve_bc_ew, ve_bc_ns,
     &                   re_bc_ns, re_bc_ew,
     &                   nx_n,ny_n,dx_n,dy_n, density, alpha)

c calculates coeficient alpha and multiplies NS BC's by it,
c so that KE density is the same accross the boundary
c Notes: add_ke2 ( in Bits and Pieces folder)

      implicit double precision (a-h,o-z)


      double precision u_bc_ew(ny_n,2), u_bc_ns(nx_n,2)
      double precision v_bc_ew(ny_n,2), v_bc_ns(nx_n,2)

      double precision re_bc_ew(ny_n,2), re_bc_ns(nx_n,2)
      double precision ue_bc_ew(ny_n,2), ue_bc_ns(nx_n,2)
      double precision ve_bc_ew(ny_n,2), ve_bc_ns(nx_n,2)


      c1 = 0.d0
      c2 = 0.d0

c length of the NS boundary
c      gamma_length = 2.d0*nx_n*dx_n + 2.d0*ny_n*dy_n

       
      do i = 1,nx_n


c integral c1
        c1 = c1 + dx_n * re_bc_ns(i,1) * 
     &     ( 0.*ue_bc_ns(i,1)**2 + ve_bc_ns(i,1)**2 )   !south
     &     + dx_n * re_bc_ns(i,2) *
     &     ( 0.*ue_bc_ns(i,2)**2 + ve_bc_ns(i,2)**2 )   !north

c integral c2

        c2 = c2 +  
     &    dx_n * ( 0.*u_bc_ns(i,1)**2 + v_bc_ns(i,1)**2 )   !south
     &  + dx_n * ( 0.*u_bc_ns(i,2)**2 + v_bc_ns(i,2)**2 )   !north
 
      enddo


      do j = 1,ny_n


c integral c1
        c1 = c1 + dy_n * re_bc_ew(j,1) * 
     &     ( ue_bc_ew(j,1)**2 + 0.*ve_bc_ew(j,1)**2 )   !west
     &    + dy_n * re_bc_ew(j,2) *
     &     ( ue_bc_ew(j,2)**2 + 0.*ve_bc_ew(j,2)**2 )   !east

c integral c2

        c2 = c2 +  
     &     dy_n * ( u_bc_ew(j,1)**2 + 0.*v_bc_ew(j,1)**2 )   !west
     &   + dy_n * ( u_bc_ew(j,2)**2 + 0.*v_bc_ew(j,2)**2 )   !east
 
      enddo


      alpha = sqrt(  c1 / (c2*density ))


c multiply  NS BC's by alpha 


      do i = 1,nx_n

        !u_bc_ns(i,1) = u_bc_ns(i,1) * alpha !south
        !u_bc_ns(i,2) = u_bc_ns(i,2) * alpha !north

        v_bc_ns(i,1) = v_bc_ns(i,1) * alpha !south
        v_bc_ns(i,2) = v_bc_ns(i,2) * alpha !north

      enddo



      do j = 1, ny_n

        u_bc_ew(j,1) = u_bc_ew(j,1) * alpha !west
        u_bc_ew(j,2) = u_bc_ew(j,2) * alpha !east
      
        !v_bc_ew(j,1) = v_bc_ew(j,1) * alpha  !west
        !v_bc_ew(j,2) = v_bc_ew(j,2) * alpha !east

      enddo


       write(*,*) 'c1,c2,alpha', c1, c2, alpha 
        
      return
      end



