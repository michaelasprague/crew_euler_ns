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

      subroutine generate_ix_iy(ix,iy,nb,nx,ny)

      implicit double precision(a-h,o-z)

! index arrays for simplifying periodic bc implementation
      integer ix ((-nb+1):(nx+nb))
      integer iy ((-nb+1):(ny+nb))

C make the ix,ixu,iy,iyv index arrays; this will simplify
C implementation of the periodic BCs
      il = 0
      ir = nx+1
      do i = 1, nb
        ix(il) = nx - (i-1)
        ix(ir) = i
        il = il - 1
        ir = ir + 1
      enddo
      do i = 1, nx
        ix(i) = i
      enddo
      il = 0
      iu = ny+1
      do i = 1, nb
        iy(il) = ny - (i-1)
        iy(iu) = i
        il = il - 1
        iu = iu + 1
      enddo
      do i = 1, ny
        iy(i) = i
      enddo

c      do i = (-nb+1), (nx+nb)
c         write(*,*) i, ix(i), iy(i)
c      enddo

      return
      end
