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

      subroutine outputvtk_ns(u, v, p, div, 
     &           xmin, ymin, dx, dy, 
     &           nx, ny, filename,
     &           umean,vmean,pmean,option)


c option 1: use given umean and vmean
c option 2: calc mean inside

      implicit double precision (a-h,o-z)

      double precision u (nx, ny)
      double precision v (nx, ny)
      double precision p (nx, ny)
      double precision div(nx, ny)
c      double precision error(nx, ny)

      integer nelem(4,nx*ny)
      integer option

      double precision xyz(2,(nx)*(ny))
      double precision u_vtk((nx)*(ny))
      double precision v_vtk((nx)*(ny))
      double precision p_vtk((nx)*(ny))
      double precision div_vtk((nx)*(ny))
c      double precision error_vtk((nx)*(ny))

      character filename*(*)

      nmax = (nx)*(ny) 
      nemax = (nx-1)*(ny-1)

      if (option .eq. 2) then
        do i = 1,nx
          do j =1,ny
            umean = umean + u(i,j)
            vmean = vmean + v(i,j)
            pmean = pmean + p(i,j)
          enddo
        enddo 

        umean = umean/float(nx*ny)
        vmean = vmean/float(nx*ny)
        pmean = pmean/float(nx*ny)
      endif


c      write(*,*) dx, dy, nmax, nemax, xmin, ymin

      y = ymin + 0.5*dy
      do j = 1, ny
        x = xmin + 0.5*dx
        do i = 1, nx

          ihat = i + (j-1)*(nx)

          xyz(1,ihat) = x
          xyz(2,ihat) = y 

! averaging to cell corners

          u_vtk(ihat) = u(i,j)

          v_vtk(ihat) = v(i,j)

          p_vtk(ihat) = p(i,j)

          div_vtk(ihat) = div(i,j)

c          error_vtk(ihat) = error(i,j)

! subtract mean
          u_vtk(ihat) = u_vtk(ihat) - umean
          v_vtk(ihat) = v_vtk(ihat) - vmean
          p_vtk(ihat) = p_vtk(ihat) - pmean

          x = x + dx
        enddo
        y = y + dy
      enddo

      ielem = 0
      do j = 1, ny -1
        do i = 1, nx -1

          ihat = i + (nx) * (j-1) 

          ielem = ielem + 1 

          nelem(1,ielem) = ihat
          nelem(2,ielem) = ihat + 1
          nelem(3,ielem) = ihat + (nx+1)
          nelem(4,ielem) = ihat + (nx+0)

        enddo
      enddo


      call vtk_scalar_out_ns(xyz, nelem, nemax, nmax, p_vtk,
     &       filename,'pressure')

      call vtk_scalar_add_ns(nmax,div_vtk,
     &       filename,'divergence',0) ! was 'divergence'

c      call vtk_scalar_add_ns(nmax,error_vtk,
c     &       filename,'error',0)

      call vtk_vector_add_ns(nmax,u_vtk,v_vtk,
     &       filename,'velocity',1)


      return
      end
