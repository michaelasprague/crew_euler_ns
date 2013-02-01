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

      subroutine l3_solve(x,d,e,f,b,nx,ny) 

      implicit double precision (a-h,o-z)

C primary variables
      double precision x(nx*ny) !solution
      double precision d(nx*ny-1) !main diagonal, 0th
      double precision e(nx*ny-2) !super diagonal, 1st
      double precision f(nx*ny-nx-1) !nx-th diagonal
      
      double precision b(nx*ny) !RHS

c      open(unit=48,file='l3.dat',status='unknown')


      x(1)=0.d0 
      n = nx*ny-1
      m = nx

      do i = 2,n+1
        x(i) = b(i)
      enddo

      do i = 1,n-m
        x(i+1) = x(i+1)/d(i)
        x(i+2) = x(i+2) - x(i+1)*e(i)
        x(i+m+1) = x(i+m+1) - x(i+1)*f(i)
      enddo

      do i = n-m+1,n-1
        x(i+1) = x(i+1)/d(i)
        x(i+2) = x(i+2) - x(i+1)*e(i)
      enddo

      i = n
        x(i+1) = x(i+1)/d(i)

        
c      do i = 1,n+1
c        write(48,*) x(i)
c     enddo

      return
      end
 
