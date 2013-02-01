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

      subroutine r_advect5(ru,ruy,rv,rvx,rrhox,rrhoy,
     &                    u,v,rho,p,dx,dy,ix,iy,nx,ny,nb)

      implicit double precision (a-h,o-z)

      double precision ru  (nx, ny)
      double precision ruy (nx, ny)
      double precision rv  (nx, ny)
      double precision rvx (nx, ny)
      double precision rrhox(nx, ny)
      double precision rrhoy(nx, ny)

      double precision u  (nx, ny)
      double precision v  (nx, ny)
      double precision rho(nx, ny)
      double precision p  (nx, ny)

! index arrays for simplifying periodic bc implementation
      integer ix ((-nb+1):(nx+nb))
      integer iy ((-nb+1):(nx+nb))

! flux terms
      do j = 1, ny
        do i = 1, nx

! cell centered
          ru(i,j) = 0.5d0*(u(i,j) + u(ix(i+1),j)) 
     &   * (
     &     +37.d0*(+0.5d0*( rho(i,j) + rho(ix(i+1),j)  ) * u(ix(i+1),j)
     &             +0.5d0*( rho(ix(i-1),j) + rho(i,j)  ) * u(i,j)
     &            )
     &     -8.d0*(+0.5d0*( rho(ix(i+2),j) + rho(ix(i+1),j))*u(ix(i+2),j)
     &            +0.5d0*( rho(ix(i-2),j) + rho(ix(i-1),j))*u(ix(i-1),j)
     &            )
     &     +1.d0*(+0.5d0*( rho(ix(i+3),j) + rho(ix(i+2),j))*u(ix(i+3),j)
     &            +0.5d0*( rho(ix(i-3),j) + rho(ix(i-2),j))*u(ix(i-2),j)
     &            )
     &      )  / 60.d0

          rv(i,j) = 0.5d0*(v(i,j) + v(i,iy(j+1))) 
     &    * (
     &     +37.d0*(+0.5d0*( rho(i,j) + rho(i,iy(j+1)) )*v(i,iy(j+1))
     &             +0.5d0*( rho(i,iy(j-1)) + rho(i,j) )*v(i,j)
     &           )
     &     -8.d0*(+0.5d0*( rho(i,iy(j+2)) + rho(i,iy(j+1)))*v(i,iy(j+2))
     &            +0.5d0*( rho(i,iy(j-2)) + rho(i,iy(j-1)))*v(i,iy(j-1))
     &           )
     &     +1.d0*(+0.5d0*( rho(i,iy(j+3)) + rho(i,iy(j+2)))*v(i,iy(j+3))
     &            +0.5d0*( rho(i,iy(j-3)) + rho(i,iy(j-2)))*v(i,iy(j-2))
     &           )
     &      ) / 60.d0  


! NS face centered
          ruy(i,j) = 0.5d0*(v(i,j) + v(ix(i-1),j)) 
     &   * (
     &  +37.d0*(+0.5d0*(rho(i,j)      +rho(ix(i-1),j)      )*u(i,j)
     &         +0.5d0*(rho(i,iy(j-1))+rho(ix(i-1),iy(j-1)))*u(i,iy(j-1))
     &            )
     &  -8.d0*(+0.5d0*(rho(i,iy(j+1))+rho(ix(i-1),iy(j+1)))*u(i,iy(j+1))
     &         +0.5d0*(rho(i,iy(j-2))+rho(ix(i-1),iy(j-2)))*u(i,iy(j-2))
     &           )
     &  +1.d0*(+0.5d0*(rho(i,iy(j+2))+rho(ix(i-1),iy(j+2)))*u(i,iy(j+2))
     &         +0.5d0*(rho(i,iy(j-3))+rho(ix(i-1),iy(j-3)))*u(i,iy(j-3))
     &           )
     &      ) / 60.d0  
 
! EW face centered
          rvx(i,j) = 0.5d0*(u(i,j) + u(i,iy(j-1))) 
     &   * (
     &   +37.d0*(+0.5d0*(rho(i,j)      +rho(i,iy(j-1))      )*v(i,j)
     &        +0.5d0*(rho(ix(i-1),j)+rho(ix(i-1),iy(j-1)))*v(ix(i-1),j)
     &            )
     & -8.d0*(+0.5d0*(rho(ix(i+1),j)+rho(ix(i+1),iy(j-1)))*v(ix(i+1),j)
     &      +0.5d0*(rho(ix(i-2),j)+rho(ix(i-2),iy(j-1)))*v(ix(i-2),j)
     &           )
     & +1.d0*(+0.5d0*(rho(ix(i+2),j)+rho(ix(i+2),iy(j-1)))*v(ix(i+2),j)
     &      +0.5d0*(rho(ix(i-3),j)+rho(ix(i-3),iy(j-1)))*v(ix(i-3),j)
     &           )
     &      )  / 60.d0


! EW face centered
          rrhox(i,j) = u(i,j) 
     &       * ( 
     &           + 37.d0 * (rho(i,j) + rho(ix(i-1),j)  )
     &           - 8.d0 * (rho(ix(i+1),j) + rho(ix(i-2),j))
     &           + 1.d0 * (rho(ix(i+2),j) + rho(ix(i-3),j))
     &         ) / 60.d0
 
! NS face centered 
          rrhoy(i,j) = v(i,j)
     &       * ( 
     &           + 37.d0 * (rho(i,j) + rho(i,iy(j-1))  )
     &           - 8.d0 * (rho(i,iy(j+1)) + rho(i,iy(j-2)))
     &           + 1.d0 * (rho(i,iy(j+2)) + rho(i,iy(j-3)))
     &         ) / 60.d0

        enddo
      enddo

! ---------------------------------------------------------------------
! modify to get fifth order
! ---------------------------------------------------------------------

      do j = 1, ny
        do i = 1, nx

! cell centered
          uavg = 0.5d0*abs(u(i,j) + u(ix(i+1),j))  
          ru(i,j) = ru(i,j) - uavg 
     &   * (
     &    +10.d0*(+0.5d0 * ( rho(i,j) + rho(ix(i+1),j)  ) * u(ix(i+1),j)
     &          -0.5d0 * ( rho(ix(i-1),j) + rho(i,j)  ) * u(i,j)
     &            )
     &   -5.d0*(+0.5d0 * ( rho(ix(i+2),j) + rho(ix(i+1),j))*u(ix(i+2),j)
     &       -0.5d0 * ( rho(ix(i-2),j) + rho(ix(i-1),j))*u(ix(i-1),j)
     &            )
     &   +1.d0*(+0.5d0 * ( rho(ix(i+3),j) + rho(ix(i+2),j))*u(ix(i+3),j)
     &        -0.5d0 * ( rho(ix(i-3),j) + rho(ix(i-2),j))*u(ix(i-2),j)
     &            )
     &      )  / 60.d0

          uavg = 0.5*abs(v(i,j) + v(i,iy(j+1))) 
          rv(i,j) = rv(i,j) - uavg
     &    * (
     &    +10.d0*(+ 0.5d0*( rho(i,j) + rho(i,iy(j+1)) )*v(i,iy(j+1))
     &        - 0.5d0*( rho(i,iy(j-1)) + rho(i,j) )*v(i,j)
     &          )
     &    -5.d0*(+ 0.5d0*( rho(i,iy(j+2)) + rho(i,iy(j+1)))*v(i,iy(j+2))
     &        - 0.5d0*( rho(i,iy(j-2)) + rho(i,iy(j-1)))*v(i,iy(j-1))
     &           )
     &    +1.d0*(+ 0.5d0*( rho(i,iy(j+3)) + rho(i,iy(j+2)))*v(i,iy(j+3))
     &          - 0.5d0*( rho(i,iy(j-3)) + rho(i,iy(j-2)))*v(i,iy(j-2))
     &           )
     &      )  / 60.d0


! NS face centered
          uavg = 0.5d0*abs(v(i,j) + v(ix(i-1),j))
          ruy(i,j) = ruy(i,j) - uavg
     &   * (
     &      +10.d0*(+0.5d0*(rho(i,j)      +rho(ix(i-1),j)      )*u(i,j)
     &        -0.5d0*(rho(i,iy(j-1))+rho(ix(i-1),iy(j-1)))*u(i,iy(j-1))
     &            )
     &  -5.d0*(+0.5d0*(rho(i,iy(j+1))+rho(ix(i-1),iy(j+1)))*u(i,iy(j+1))
     &      -0.5d0*(rho(i,iy(j-2))+rho(ix(i-1),iy(j-2)))*u(i,iy(j-2))
     &           )
     &  +1.d0*(+0.5d0*(rho(i,iy(j+2))+rho(ix(i-1),iy(j+2)))*u(i,iy(j+2))
     &      -0.5d0*(rho(i,iy(j-3))+rho(ix(i-1),iy(j-3)))*u(i,iy(j-3))
     &           )
     &      )  / 60.d0
 
! EW face centered
          uavg = 0.5d0*abs(u(i,j) + u(i,iy(j-1)))
          rvx(i,j) = rvx(i,j) - uavg 
     &   * (
     &   +10.d0*(+0.5d0*(rho(i,j)      +rho(i,iy(j-1))      )*v(i,j)
     &        -0.5d0*(rho(ix(i-1),j)+rho(ix(i-1),iy(j-1)))*v(ix(i-1),j)
     &            )
     & -5.d0*(+0.5d0*(rho(ix(i+1),j)+rho(ix(i+1),iy(j-1)))*v(ix(i+1),j)
     &      -0.5d0*(rho(ix(i-2),j)+rho(ix(i-2),iy(j-1)))*v(ix(i-2),j)
     &           )
     & +1.d0*(+0.5d0*(rho(ix(i+2),j)+rho(ix(i+2),iy(j-1)))*v(ix(i+2),j)
     &      -0.5d0*(rho(ix(i-3),j)+rho(ix(i-3),iy(j-1)))*v(ix(i-3),j)
     &           )
     &      )  / 60.d0


! EW face centered
          uavg = abs(u(i,j))
          rrhox(i,j) = rrhox(i,j) - uavg 
     &       * ( 
     &           + 10.d0 * (rho(i,j) - rho(ix(i-1),j)  )
     &           - 5.d0 * (rho(ix(i+1),j) - rho(ix(i-2),j))
     &           + 1.d0 * (rho(ix(i+2),j) - rho(ix(i-3),j))
     &         ) / 60.d0
 
! NS face centered 
          uavg = abs(v(i,j))
          rrhoy(i,j) = rrhoy(i,j) - uavg
     &       * ( 
     &           + 10.d0 * (rho(i,j) - rho(i,iy(j-1))  )
     &           - 5.d0 * (rho(i,iy(j+1)) - rho(i,iy(j-2)))
     &           + 1.d0 * (rho(i,iy(j+2)) - rho(i,iy(j-3)))
     &         ) / 60.d0

        enddo
      enddo


      end
