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

      subroutine init_euler(u,v,rho,x,xc,y,yc,ix,iy,nx,ny,nb,
     &                      vortex_center, umean, vmean, rhomean)

      implicit double precision (a-h,o-z)

      double precision u  (nx, ny)
      double precision v  (nx, ny)
      double precision rho(nx, ny)

      double precision x (nx)
      double precision xc(nx)
      double precision y (ny)
      double precision yc(ny)

      double precision vortex_center(2)


! index arrays for simplifying periodic bc implementation
      integer ix ((-nb+1):(nx+nb))
      integer iy ((-nb+1):(ny+nb))

! initial condition taken from 
!  Y.C. Zhou and G.W. Wei, High resolution conjugate filters for
!  the simulations of flows, JCP 189, 2003, 159-179

      gam = 1.4d0  

      dlambda = 5.d0 ! strength of the vortex

      eta = 1.d0

      pi = acos(-1.d0)

      x0 = vortex_center(1) !-5.d0
      y0 = vortex_center(2) !-5.d0

c      umean = 1.d0
c      vmean = 1.d0
c      rhomean = 1.d0

      do j = 1, ny
        do i = 1, nx

          ru2 = (x(i)-x0)**2 + (yc(j)-y0)**2 

          rv2 = (xc(i)-x0)**2 + (y(j)-y0)**2 

          rc2 = (xc(i)-x0)**2 + (yc(j)-y0)**2 

          u  (i,j) = umean - (0.5d0*dlambda / pi) 
     &             * (yc(j) - y0)*exp( eta * (1.d0-ru2))

          v  (i,j) = vmean + (0.5d0*dlambda / pi) 
     &             * (xc(i) - x0)*exp( eta * (1.d0-rv2))

          rho(i,j) = (rhomean - (gam-1.d0) * dlambda**2 
     &                    * exp( 2.d0 * eta * (1.d0-rc2))
     &                / (16.d0 * eta * gam * pi**2))**(1.d0/(gam-1.d0))

        enddo
      enddo

      return
      end
