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

      subroutine update_euler(t, u, v, r, p,
     &                        rux, ruy, rvy, rvx, rrx, rry,
     &                        dtau,iacoustic, dx, dy, u_max,
     &                        ix, iy, nb, nx, ny, k)

      implicit double precision(a-h,o-z)

! solution at time t 
      double precision u (nx, ny)
      double precision v (nx, ny)
      double precision r (nx, ny)
      double precision p (nx, ny)

! flux terms
      double precision rux (nx, ny)
      double precision rvy (nx, ny)
      double precision ruy (nx, ny)
      double precision rvx (nx, ny)
      double precision rrx (nx, ny)
      double precision rry (nx, ny)

! solution at time t + tau (acoustic solution)
      double precision ua (nx, ny)
      double precision va (nx, ny)
      double precision ra (nx, ny)
      double precision pa (nx, ny)

! solution at time t + tau + dtau(acoustic solution)
      double precision ua_new (nx, ny)
      double precision va_new (nx, ny)
      double precision ra_new (nx, ny)
      double precision pa_new (nx, ny)

! index arrays for simplifying periodic bc implementation
      integer ix ((-nb+1):(nx+nb))
      integer iy ((-nb+1):(ny+nb))

      integer iacoustic(3)
      double precision  dtau(3)


      if (k.eq.1) then
        call r_advect5(rux,ruy,rvy,rvx,rrx,rry,
     &                 u,v,r,p,dx,dy,
     &                 ix,iy,nx,ny,nb)
      endif


      do krk = 1, 3  ! rk3 substepping
        
        call copy(u, ua, nx, ny)
        call copy(v, va, nx, ny)
        call copy(r, ra, nx, ny)
        call copy(p, pa, nx, ny)

 
        tau = t ! acoustic time -- maybe this should be zero ??
        do kac = 1, iacoustic(krk) ! acoustic substepping
           
          do j = 1, ny 
            do i = 1, nx 

c----------------------------------------------------
c         start the business
c----------------------------------------------------

             ua_new(i,j) = 
     &           + ua(i,j)*0.5d0 * ( ra(i,j) + ra(ix(i-1),j) )
     &           - (dtau(krk)/dx) * (pa(i,j) - pa(ix(i-1),j)  )
     &           - (dtau(krk)/dx) * (rux(i,j) - rux(ix(i-1),j))
     &           - (dtau(krk)/dy) * (ruy(i,iy(j+1)) - ruy(i,j))
c artifical divergence damping
c    &- 0.025*dtau(krk)*ra(i,j)*(ua(ix(i+1),j) - ua(ix(i-1),j))/dx**2

             va_new(i,j) = 
     &           + va(i,j)*0.5d0*( ra(i,j) + ra(i,iy(j-1)) )
     &           - (dtau(krk)/dy) * (pa(i,j)  - pa(i,iy(j-1)) )
     &           - (dtau(krk)/dy) * (rvy(i,j) - rvy(i,iy(j-1)))
     &           - (dtau(krk)/dx) * (rvx(ix(i+1),j) - rvx(i,j))
c    &- 0.025*dtau(krk)*ra(i,j)*(va(i,iy(j+1)) - va(i,iy(j-1)))/dy**2

               ra_new(i,j) = ra(i,j) 
     &           - (dtau(krk)/dx) * (rrx(ix(i+1),j) - rrx(i,j))
     &           - (dtau(krk)/dy) * (rry(i,iy(j+1)) - rry(i,j))



c----------------------------------------------------
c         finish the business
c----------------------------------------------------

            enddo 
          enddo 

          call copy(ua_new, ua, nx, ny)
          call copy(va_new, va, nx, ny)
          call copy(ra_new, ra, nx, ny)

          call deflux(ua, va, ra, ix, iy, nx, ny, nb)
 
          call update_pressure(pa, ra, t, nx, ny)

!         call p_acoustic(pa,p,ra,r,nx,ny)
 
          tau = tau + dtau(krk)
 
        enddo


!  here we have the acoustic values at the end of substeps

        call r_advect5(rux, ruy, rvy, rvx, rrx, rry,
     &                 ua, va, ra, pa, dx, dy,
     &                 ix, iy, nx, ny, nb)

       enddo

       call copy(ua, u, nx, ny)
       call copy(va, v, nx, ny)
       call copy(ra, r, nx, ny)
       call copy(pa, p, nx, ny)

! check to ensure u is less than u_max specified above; trying to understand
! instabilities that are being exhibited; 12/9/2011a
       do j = 1,ny
         do i = 1,nx
           if (u(i,j) .gt. u_max) write(*,*) 'u = ',u(i,j)
           if (v(i,j) .gt. u_max) write(*,*) 'v = ',v(i,j)
         enddo
       enddo

       return
       end

