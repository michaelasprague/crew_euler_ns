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

      subroutine beta(w, ew, ns, vn, An, Hnvs, Hnvn, 
     &                  bn, bn1, fs, fn, 
     &                  nx,ny,dx,dy,dt,density, k)


c calculates quantitiesat the n-s edges and e-w edges
c see pg.12 for notes
c i.e. calculates terms on the RHS of press-Laplace eq.



      implicit double precision (a-h, o-z)

      double precision ew(nx-1,ny)
      double precision ns(nx,ny-1)
      double precision vn(nx,ny) ! u_star 
  
      double precision Hnvs(nx,ny) ! H^n(v^*)
      double precision Hnvn(nx,ny) ! H^n(v^n)
      double precision An(nx,ny) ! A^* times ones 

      double precision bn (2*nx + 2*ny - 4, 2)
      double precision bn1(2*nx + 2*ny - 4, 2)

      double precision fs(nx,ny,2), fn(nx,ny,2) ! forcing at u_star and u_n

      double precision w(nx,ny) ! inside array


      a = (dx*dy)/dt
 
 
c calculating "two" at the centers

      do j = 1,ny
        do i = 1,nx

          w(i,j) =  
     &            (a + 0.5d0*An(i,j)) * vn(i,j)
     &            + 0.5d0 * ( Hnvs(i,j) + Hnvn(i,j) )
     &            + 0.5d0 * dx * dy * ( fs(i,j,k) + fn(i,j,k) ) 

        enddo
      enddo

c South
      do i = 1,nx
        w(i,1) = w(i,1) + 0.5d0 * ( bn(i,k) + bn1(i,k) )
      enddo


c East
      do i = nx+1, nx+ny-1
         j = i - nx + 1 ! so j = 2,ny
         w(nx,j) = w(nx,j) + 0.5d0 * ( bn(i,k) + bn1(i,k) )
      enddo

c North 
      do i = nx+ny, 2*nx+ny-2
         j = -(i-2*nx-ny+1) ! so j = nx-1,1
         w(j,ny) = w(j,ny) + 0.5d0 * ( bn(i,k) + bn1(i,k) )
      enddo

c West
      do i = 2*nx+ny-1, 2*nx+2*ny-4
        j = -(i-2*nx-2*ny+2) ! so j = ny-1, 2
        w(1,j) = w(1,j) + 0.5d0 * ( bn(i,k) + bn1(i,k) )
      enddo

      do j = 1,ny
        do i = 1, nx
          w(i,j) = w(i,j) / ( a - 0.5d0*An(i,j) )
        enddo
      enddo

c interpolating 

      if (k.eq.1) then
        do j = 1,ny
          do i = 1, nx-1

            ew(i,j) = 0.5d0* dy * ( w(i,j) + w(i+1,j) )

          enddo
        enddo
      else
        do j = 1,ny-1
          do i = 1, nx

            ns(i,j) = 0.5d0 * dx * ( w(i,j) + w(i,j+1) )

          enddo
        enddo
      endif

      return
      end
