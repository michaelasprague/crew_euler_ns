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

      subroutine chol(d,e,f,
     &                A_ew,A_ns,dx,dy,dt,nx,ny,density) 

c Incomplete Cholesky preconditioner of press matrix.
c Obtains diagonals of press matrix using p_mat_vec
c (with first row and columd removed i.e. fixed p(1,1) =0 )
c Then performs incomplete Chol decomposition.



      implicit double precision (a-h,o-z)

      double precision A_ew((nx-1)*ny)
      double precision A_ns(nx*(ny-1))
C chol returns: 
      double precision d(nx*ny-1) !main diagonal, 0th
      double precision e(nx*ny-2) !super diagonal, 1st
      double precision f(nx*ny-nx-1) !nx-th diagonal
      
c      double precision one(nx*ny) !basis vector 
c      double precision col(nx*ny) !column of P 
     
c      external p_mat_vec
c      external initialize 
c      external get_p_diag

c      open(unit=45,file='d.dat',status='unknown')
c      open(unit=46,file='e.dat',status='unknown')
c      open(unit=47,file='f.dat',status='unknown')
      n = nx*ny
      m = nx
 

cc-----------getting diagonals of P (smart way)-----------

      call get_p_diag(d,e,f,A_ew,A_ns,
     $                nx,ny,dx,dy,dt,density)
 
c-------------cholesky incomplete ------------------------
c       write(*,*) 4.d0**0.5d0
      n = n-1

      do k = 1,n-m
        d(k) = d(k)**0.5d0
        e(k) = e(k)/d(k)
        f(k) = f(k)/d(k)
   
        d(k+1) = d(k+1) - e(k)*e(k)
        d(k+m) = d(k+m) - f(k)*f(k)
      enddo
 
      do k = n-m+1,n-1
        d(k) = d(k)**0.5d0
        e(k) = e(k)/d(k)
        d(k+1) = d(k+1) - e(k)*e(k)
      enddo

      k=n
      d(k) = d(k)**0.5d0


      
c      do i = 1,n-1
c         write(45,*) d(i)
c      enddo
c
c      do i = 1,n-2
c         write(46,*) e(i)
c      enddo
c
c      do i = 1,n-m-1
c        write(47,*) f(i)
c      enddo

      return
      end

 
