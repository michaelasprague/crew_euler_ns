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

      subroutine outputvtk(u,v,rho,p,xmin,ymin,
     &           dx,dy,ix,iy,nx,ny,nb,filename,option)

      implicit double precision (a-h,o-z)

      double precision u  (nx, ny)
      double precision v  (nx, ny)
      double precision rho(nx, ny)
      double precision p(nx, ny)
    


! index arrays for simplifying periodic bc implementation
      integer ix ((-nb+1):(nx+nb))
      integer iy ((-nb+1):(ny+nb))
      integer option ! 1 for mean = 1, 2 for mean = 0


      integer nelem(4,nx*ny)
      double precision xyz(2,(nx+1)*(ny+1))
      double precision u_vtk((nx+1)*(ny+1))
      double precision v_vtk((nx+1)*(ny+1))
      double precision rho_vtk((nx+1)*(ny+1))
      double precision p_vtk((nx+1)*(ny+1))


      character filename*(*)
      logical mass
c if (mass) then output mass flux perturbation, as opposed to, velocity
c perturbation


      if (option .eq. 1) then
         umean = 1.d0
         vmean = 1.d0
         rhomean = 1.d0
         pmean = 1.d0
      elseif (option .eq. 2) then
         umean = 0.d0
         vmean = 0.d0
         rhomean = 0.d0
         pmean = 0.d0
      endif

      nmax = (nx+1)*(ny+1) 
      nemax = nx*ny

      y = ymin
      do j = 1, ny+1
        x = xmin
        do i = 1, nx+1
          ihat = i + (j-1)*(nx+1)
          xyz(1,ihat) = x
          xyz(2,ihat) = y 

! averaging to cell corners

          u_vtk(ihat) = 0.5*(u(ix(i),iy(j-1)) + u(ix(i),iy(j)))

          v_vtk(ihat) = 0.5*(v(ix(i),iy(j)) + v(ix(i-1),iy(j)))

          rho_vtk(ihat) = 0.25*(
     &                    rho(ix(i),iy(j)) + 
     &                    rho(ix(i-1),iy(j)) + 
     &                    rho(ix(i-1),iy(j-1)) + 
     &                    rho(ix(i),iy(j-1)) )

          p_vtk(ihat) = 0.25*(
     &                    p(ix(i),iy(j)) + 
     &                    p(ix(i-1),iy(j)) + 
     &                    p(ix(i-1),iy(j-1)) + 
     &                    p(ix(i),iy(j-1)) )


! subtract mean
           
          u_vtk(ihat) = u_vtk(ihat) - umean
          v_vtk(ihat) = v_vtk(ihat) - vmean

          rho_vtk(ihat) = rho_vtk(ihat) - rhomean
          p_vtk(ihat) = p_vtk(ihat) - pmean


          x = x + dx
        enddo
        y = y + dy
      enddo

      ielem = 0
      do j = 1, ny
        do i = 1, nx

          ihat = i + (nx+1) * (j-1) 

          ielem = ielem + 1 

          nelem(1,ielem) = ihat
          nelem(2,ielem) = ihat + 1
          nelem(3,ielem) = ihat + (nx+2)
          nelem(4,ielem) = ihat + (nx+1)

        enddo
      enddo


      call vtk_scalar_out(xyz,nelem,nemax,nmax,rho_vtk,
     &       filename,'density')

      call vtk_scalar_add(xyz,nelem,nemax,nmax,p_vtk,
     &         filename,'pressure',0)

      call vtk_vector_add(xyz,nelem,nemax,nmax,u_vtk,v_vtk,
     &         filename,'velocity',1) !1



      return
      end
