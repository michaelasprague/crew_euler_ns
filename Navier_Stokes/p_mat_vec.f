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

      subroutine p_mat_vec(w,p,A_ew,A_ns,A_e,A_n,
     &                     nx,ny,dx,dy,dt,density,neum)

      implicit double precision (a-h,o-z)

      double precision p(nx,ny), w(nx,ny)
c     double precision As(nx,ny)
      double precision A_ew(nx-1,ny)
      double precision A_ns(nx,ny-1)
      double precision A_e(ny) !(a_tilde)^-1 at centers along East edge
      double precision A_n(nx) !(a_tilde)^-1 at centers along North edge
 
      double precision d(nx*ny-1) !main diagonal, 0th
      double precision e(nx*ny-2) !super diagonal, 1st
      double precision f(nx*ny-nx-1) !ny-th diagonal

      logical neum
      logical symetric
     

      if (neum) symetric = .true.! .false.
      if (.not. neum) symetric = .true.

      

      dx2 = dx*dx
      dy2 = dy*dy
      a = (dx*dy)/dt
  
c Inside
!$omp parallel default(shared) private(i,j) 
!$omp do schedule(dynamic)
      do j = 2, ny-1
        do i = 2, nx-1
           w(i,j) = -( (A_ew(i,j)+A_ew(i-1,j))*dy2
     &                +(A_ns(i,j)+A_ns(i,j-1))*dx2 ) * p(i,j)
     &            + dy2*A_ew(i,j)*p(i+1,j) + dy2*A_ew(i-1,j)*p(i-1,j) 
     &            + dx2*A_ns(i,j)*p(i,j+1) + dx2*A_ns(i,j-1)*p(i,j-1)

        enddo
      enddo
!$omp end parallel

c South
      j = 1
      do i = 2, nx-1
        w(i,j) = -( (A_ew(i,j) + A_ew(i-1,j))*dy2
     &             +(A_ns(i,j) + 0.d0       )*dx2 ) * p(i,j)
     &            + dy2*A_ew(i,j) * p(i+1,j) + dy2*A_ew(i-1,j)*p(i-1,j) 
     &            + dx2*A_ns(i,j) * p(i,j+1) + 0.d0 

      enddo

c North
      j = ny 
      do i = 2, nx-1

        if (.not. symetric) then

        w(i,j) = -( (A_ew(i,j)+A_ew(i-1,j))*dy2
     &             +(-0.5d0*A_n(i) +A_ns(i,j-1))*dx2 ) * p(i,j)
     &            + dy2*A_ew(i,j)*p(i+1,j) + dy2*A_ew(i-1,j)*p(i-1,j)
     &            + 0.d0    + dx2*(A_ns(i,j-1)-0.5d0*A_n(i))*p(i,j-1)

        else

        w(i,j) = -( (A_ew(i,j)+A_ew(i-1,j))*dy2
     &             +(0.d0     +A_ns(i,j-1))*dx2 ) * p(i,j)
     &            + dy2*A_ew(i,j)*p(i+1,j) + dy2*A_ew(i-1,j)*p(i-1,j) 
     &            + 0.d0                   + dx2*A_ns(i,j-1)*p(i,j-1)

        endif

      enddo

c West
      i = 1
      do j = 2, ny-1
        w(i,j) = -( (A_ew(i,j)+0.d0       )*dy2
     &             +(A_ns(i,j)+A_ns(i,j-1))*dx2 ) * p(i,j)
     &            + dy2*A_ew(i,j)*p(i+1,j) + 0.d0 
     &            + dx2*A_ns(i,j)*p(i,j+1) + dx2*A_ns(i,j-1)*p(i,j-1)
      enddo

c East

      i = nx
      do j = 2, ny-1

        if (.not. symetric) then

        w(i,j) = -( (-0.5d0*A_e(j)+A_ew(i-1,j))*dy2
     &             +(A_ns(i,j)+A_ns(i,j-1))*dx2 ) * p(i,j)
     &            + 0.d0    + dy2*(A_ew(i-1,j)-0.5d0*A_e(j))*p(i-1,j)
     &            + dx2*A_ns(i,j)*p(i,j+1) + dx2*A_ns(i,j-1)*p(i,j-1)
        
        else

        w(i,j) = -( (0.d0     +A_ew(i-1,j))*dy2
     &             +(A_ns(i,j)+A_ns(i,j-1))*dx2 ) * p(i,j)
     &            + 0.d0                   + dy2*A_ew(i-1,j)*p(i-1,j) 
     &            + dx2*A_ns(i,j)*p(i,j+1) + dx2*A_ns(i,j-1)*p(i,j-1)

        endif

      enddo


c Corners:

c SW

      i = 1
      j = 1

       w(i,j) = -( (A_ew(i,j)+0.d0       )*dy2
     &            +(A_ns(i,j)+0.d0       )*dx2 ) * p(i,j)
     &            + dy2*A_ew(i,j)*p(i+1,j) + 0.d0 
     &            + dx2*A_ns(i,j)*p(i,j+1) + 0.d0 

c SE

      i = nx 
      j = 1

       if (.not. symetric) then

       w(i,j) = -( (-0.5d0*A_e(j)     +A_ew(i-1,j))*dy2
     &            +(A_ns(i,j)+0.d0        )*dx2 ) * p(i,j)
     &            + 0.d0   + dy2*(A_ew(i-1,j)-0.5d0*A_e(j))*p(i-1,j) 
     &            + dx2*A_ns(i,j)*p(i,j+1) + 0.d0 

       else

       w(i,j) = -( (0.d0       +A_ew(i-1,j))*dy2
     &            +(A_ns(i,j)+0.d0        )*dx2 ) * p(i,j)
     &            + 0.d0        + dy2*A_ew(i-1,j)*p(i-1,j)
     &            + dx2*A_ns(i,j)*p(i,j+1) + 0.d0

       endif


c NW

      i = 1
      j = ny

       if (.not. symetric) then
    
       w(i,j) = -( (A_ew(i,j)+0.d0       )*dy2
     &            +(-0.5d0*A_n(i)   +A_ns(i,j-1))*dx2 ) * p(i,j)
     &            + dy2*A_ew(i,j)*p(i+1,j) + 0.d0
     &            + 0.d0     + dx2*(A_ns(i,j-1)-0.5d0*A_n(i))*p(i,j-1)

       else

       w(i,j) = -( (A_ew(i,j)+0.d0       )*dy2
     &            +(0.d0     +A_ns(i,j-1))*dx2 ) * p(i,j)
     &            + dy2*A_ew(i,j)*p(i+1,j) + 0.d0 
     &            + 0.d0                   + dx2*A_ns(i,j-1)*p(i,j-1)

       endif

c NE

      i = nx
      j = ny

       if (.not. symetric) then

       w(i,j) = -( (-0.5d0*A_e(i)     +A_ew(i-1,j))*dy2
     &            +(-0.5d0*A_n(j)     +A_ns(i,j-1))*dx2 ) * p(i,j)
     &            + 0.d0     + dy2*(A_ew(i-1,j)-0.5d0*A_e(i))*p(i-1,j)
     &            + 0.d0     + dx2*(A_ns(i,j-1)-0.5d0*A_n(j))*p(i,j-1)

       else

       w(i,j) = -( (0.d0     +A_ew(i-1,j))*dy2
     &            +(0.d0     +A_ns(i,j-1))*dx2 ) * p(i,j)
     &            + 0.d0                   + dy2*A_ew(i-1,j)*p(i-1,j) 
     &            + 0.d0                   + dx2*A_ns(i,j-1)*p(i,j-1)

       endif

c 
      do j = 1, ny
        do i = 1, nx
          w(i,j) = w(i,j) / (-1.d0 * density) !negative for pos def matrix
        enddo                                 !rhs is also multiplied by -1
      enddo

      return
      end
