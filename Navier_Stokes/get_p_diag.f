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

      subroutine get_p_diag(d,e,f,A_ew,A_ns,
     &                     nx,ny,dx,dy,dt,density)

c get the diagonals of P matrix (with first col and row removed)
c d - main, e-super, f- nx-th 

      implicit double precision (a-h,o-z)

c      double precision p(nx,ny), w(nx,ny)
c     double precision As(nx,ny)
      double precision A_ew(nx-1,ny)
      double precision A_ns(nx,ny-1)

      double precision d(nx*ny-1) !main diagonal, 0th
      double precision e(nx*ny-2) !super diagonal, 1st
      double precision f(nx*ny-nx-1) !ny-th diagonal


      dx2 = dx*dx
      dy2 = dy*dy
      a = (dx*dy)/dt

c Inside
      do j = 2, ny-1
        do i = 2, nx-1

           d(i+(j-1)*nx-1) =  (A_ew(i,j)+A_ew(i-1,j))*dy2
     &                       +(A_ns(i,j)+A_ns(i,j-1))*dx2  

           e(i+(j-1)*nx-1-1) = -dy2*A_ew(i-1,j)
           f(i+(j-1)*nx-1-nx) = -dx2*A_ns(i,j-1)


c           w(i,j) = -( (A_ew(i,j)+A_ew(i-1,j))*dy2
c     &                +(A_ns(i,j)+A_ns(i,j-1))*dx2 ) * p(i,j)
c     &            + dy2*A_ew(i,j)*p(i+1,j) + dy2*A_ew(i-1,j)*p(i-1,j) 
c     &            + dx2*A_ns(i,j)*p(i,j+1) + dx2*A_ns(i,j-1)*p(i,j-1)

        enddo
      enddo

c South
      j = 1
      do i = 3, nx-1

        d(i+(j-1)*nx-1) = (A_ew(i,j) + A_ew(i-1,j))*dy2
     &                   +(A_ns(i,j) + 0.d0       )*dx2

        e(i+(j-1)*nx-1-1) = -dy2*A_ew(i-1,j) 
c        f(i+(j-1)*nx-1) = 0.d0 

c        w(i,j) = -( (A_ew(i,j) + A_ew(i-1,j))*dy2
c     &             +(A_ns(i,j) + 0.d0       )*dx2 ) * p(i,j)
c     &            + dy2*A_ew(i,j) * p(i+1,j) + dy2*A_ew(i-1,j)*p(i-1,j) 
c     &            + dx2*A_ns(i,j) * p(i,j+1) + 0.d0 

      enddo
      i=2
        d(i+(j-1)*nx-1) = (A_ew(i,j) + A_ew(i-1,j))*dy2
     &                   +(A_ns(i,j) + 0.d0       )*dx2

c        e(i+(j-1)*nx-1-1) = -dy2*A_ew(i-1,j)



c North
      j = ny 
      do i = 2, nx-1
      
        d(i+(j-1)*nx-1) = (A_ew(i,j)+A_ew(i-1,j))*dy2
     &                   +(0.d0     +A_ns(i,j-1))*dx2
     
        e(i+(j-1)*nx-1-1) = - dy2*A_ew(i-1,j)

        f(i+(j-1)*nx-1-nx) = -A_ns(i,j-1)*dx2
  

c        w(i,j) = -( (A_ew(i,j)+A_ew(i-1,j))*dy2
c     &             +(0.d0     +A_ns(i,j-1))*dx2 ) * p(i,j)
c     &            + dy2*A_ew(i,j)*p(i+1,j) + dy2*A_ew(i-1,j)*p(i-1,j) 
c     &            + 0.d0                   + dx2*A_ns(i,j-1)*p(i,j-1)
      enddo

c West
      i = 1
      do j = 2, ny-1

        d(i+(j-1)*nx-1) = (A_ew(i,j)+0.d0       )*dy2
     &                   +(A_ns(i,j)+A_ns(i,j-1))*dx2
        e(i+(j-1)*nx-1-1) = 0.d0    
        f(i+(j-1)*nx-1-nx) = - dx2*A_ns(i,j-1)
 
c        w(i,j) = -( (A_ew(i,j)+0.d0       )*dy2
c     &             +(A_ns(i,j)+A_ns(i,j-1))*dx2 ) * p(i,j)
c     &            + dy2*A_ew(i,j)*p(i+1,j) + 0.d0 
c     &            + dx2*A_ns(i,j)*p(i,j+1) + dx2*A_ns(i,j-1)*p(i,j-1)
      enddo

c East

      i = nx
      do j = 2, ny-1
        
        d(i+(j-1)*nx-1) = (0.d0     +A_ew(i-1,j))*dy2
     &             +(A_ns(i,j)+A_ns(i,j-1))*dx2 
        e(i+(j-1)*nx-1-1) = -dy2*A_ew(i-1,j)
        f(i+(j-1)*nx-1-nx) = -dx2*A_ns(i,j-1)


c        w(i,j) = -( (0.d0     +A_ew(i-1,j))*dy2
c     &             +(A_ns(i,j)+A_ns(i,j-1))*dx2 ) * p(i,j)
c     &            + 0.d0                   + dy2*A_ew(i-1,j)*p(i-1,j) 
c     &            + dx2*A_ns(i,j)*p(i,j+1) + dx2*A_ns(i,j-1)*p(i,j-1)

      enddo


c Corners:

c SW

c      i = 1
c      j = 1
c
c       w(i,j) = -( (A_ew(i,j)+0.d0       )*dy2
c     &            +(A_ns(i,j)+0.d0       )*dx2 ) * p(i,j)
c     &            + dy2*A_ew(i,j)*p(i+1,j) + 0.d0 
c     &            + dx2*A_ns(i,j)*p(i,j+1) + 0.d0 

c SE

      i = nx 
      j = 1

       d(i+(j-1)*nx-1) = (0.d0     +A_ew(i-1,j))*dy2
     &            +(A_ns(i,j)+0.d0        )*dx2    

       e(i+(j-1)*nx-1-1) = - dy2*A_ew(i-1,j)
c       f(i+(j-1)*nx-1) = 0.d0

c       w(i,j) = -( (0.d0     +A_ew(i-1,j))*dy2
c     &            +(A_ns(i,j)+0.d0        )*dx2 ) * p(i,j)
c     &            + 0.d0                   + dy2*A_ew(i-1,j)*p(i-1,j) 
c     &            + dx2*A_ns(i,j)*p(i,j+1) + 0.d0 

c NW

      i = 1
      j = ny

       d(i+(j-1)*nx-1) = (A_ew(i,j)+0.d0       )*dy2
     &            +(0.d0     +A_ns(i,j-1))*dx2 

       e(i+(j-1)*nx-1-1) = 0.d0
       f(i+(j-1)*nx-1-nx) = - dx2*A_ns(i,j-1)
 
c       w(i,j) = -( (A_ew(i,j)+0.d0       )*dy2
c     &            +(0.d0     +A_ns(i,j-1))*dx2 ) * p(i,j)
c     &            + dy2*A_ew(i,j)*p(i+1,j) + 0.d0 
c     &            + 0.d0                   + dx2*A_ns(i,j-1)*p(i,j-1)

c NE

      i = nx
      j = ny
       
       d(i+(j-1)*nx-1) = (0.d0     +A_ew(i-1,j))*dy2
     &            +(0.d0     +A_ns(i,j-1))*dx2 

       e(i+(j-1)*nx-1-1) = -dy2*A_ew(i-1,j)
       f(i+(j-1)*nx-1-nx) = -dx2*A_ns(i,j-1)

c       w(i,j) = -( (0.d0     +A_ew(i-1,j))*dy2
c     &            +(0.d0     +A_ns(i,j-1))*dx2 ) * p(i,j)
c     &            + 0.d0                   + dy2*A_ew(i-1,j)*p(i-1,j) 
c     &            + 0.d0                   + dx2*A_ns(i,j-1)*p(i,j-1)

c      do j = 1, ny
c        do i = 1, nx
c          w(i,j) = w(i,j) / (-1.d0 * density) !negative for pos def
c        enddo
c      enddo

      return
      end
