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

      subroutine bicgstab(u, rhs,
     &                    ew_rc, ns_rc,
     &                    viscosity,
     &                    dx,dy,dt,nx,ny,dmat_vec,ierr,kfinal,tol,
     &                    neum,conv)

      implicit double precision (a-h,o-z)

C primary variables
c      double precision u(nx*ny), v(nx*ny), p(nx*ny)
      double precision u(nx*ny)

      double precision ew_rc(nx+1,ny), ns_rc(nx,ny+1)

C RHS forcing in pressure equation
      double precision rhs(nx*ny)

C   Conjugate gradient arrays
      double precision w(nx*ny), r(nx*ny), z(nx*ny), p(nx*ny), s(nx*ny)
      double precision ws(nx*ny), rhat(nx*ny)

      external dmat_vec
 
      logical neum
      logical  conv

      nxny = nx*ny

      kmax = nxny*100

      ierr = 0

c      tol = 1.d-5
        
      eps = 1.d-30

      call dmat_vec(w, u,
     &              ew_rc, ns_rc,
     &              nx, ny, dx, dy, dt, viscosity, neum,conv)

      rmax = 0.d0
      rhs_max = 0.d0
      r2norm = 0.d0
      rhs2norm = 0.d0

      do i = 1,nxny

        r(i) = rhs(i) - w(i)

c        write(*,*) i,rhs(i),w(i),r(i)

        p(i) = r(i)

        rhat(i) = r(i)

        if (abs(r(i)) .gt. rmax)  rmax = abs(r(i))
        if (abs(rhs(i)) .gt. rhs_max) rhs_max = abs(rhs(i))
        r2norm = r2norm + r(i)**2.d0
        rhs2norm = rhs2norm + rhs(i)**2.d0
      enddo

      r2norm = r2norm**0.5d0
      rhs2norm = rhs2norm**0.5d0

      rhs_max = 1.d0 ! i.e. no rel tolerance cond
      
c      write(*,*) 'RHS_MAX', rhs_max

c - this makes me nervous -- mas 10 jan 2012
      if (rmax .lt. tol*rhs_max) then
c      if (r2norm .lt. tol*rhs2norm) then
        k = 0
        do i = 1, nxny
          p(i) = u(i)
        enddo
        goto 5   
      endif

      alphatop = 0.d0

      do i = 1,nxny
        alphatop = alphatop + rhat(i) * r(i)
      enddo

! check to ensure inner product rhat r is not zero
      if (abs(alphatop) .lt. eps) then
         write(*,*) 'rmax = ',rmax
         write(*,*) 'alphatop = ',alphatop
         stop 'inner product zero; choose different rhat'
      endif 

      do k = 1,kmax

        call dmat_vec(w,p,
     &               ew_rc,ns_rc,
     &               nx,ny,dx,dy,dt,viscosity, neum,conv)

        alphabot = 0.d0

        do i = 1,nxny
          alphabot = alphabot + rhat(i) * w(i)
        enddo

        !if (abs(alphabot) .lt. eps .and. alphatop .gt. eps) then
        if (abs(alphabot) .lt. eps ) then
          write(*,*) 'alphabot  is zero'
          write(*,*) 'alphatop  is', alphatop
          stop
        endif

        alpha_cg = alphatop / alphabot
 
        do i = 1, nxny
          s(i) = r(i) - alpha_cg * w(i)
        enddo

        call dmat_vec(ws,s,
     &               ew_rc,ns_rc,
     &               nx,ny,dx,dy,dt,viscosity, neum,conv)

        omegatop = 0.d0
        omegabot = 0.d0
        do i = 1, nxny
          omegatop = omegatop + ws(i)*s(i)
          omegabot = omegabot + ws(i)*ws(i)
        enddo

        if (abs(omegabot) .lt. eps) then
          omega = 0.d0
          stop 'omegabot too small'
        else
          omega = omegatop / omegabot
        endif

        rmax = 0.d0
        r2norm = 0.d0
 
        do i = 1, nxny

          u(i) = u(i) + alpha_cg*p(i) + omega * s(i)

          r(i) = s(i) - omega*ws(i)  ! r^(k+1) = r^(k) - alpha A p^(k)

          if ( abs(r(i)) .gt. rmax)  rmax = abs(r(i))
          r2norm = r2norm + r(i)**2.d0
        enddo
          r2norm = r2norm**0.5d0

        if (rmax .lt. tol*rhs_max .and. k .ge. 1)  goto 5 ! converged
c        if (r2norm .lt. tol*rhs2norm .and. k .ge. 1)  goto 5 ! converged

        if (rmax .gt. 1d6) then
          write(*,*) '---------- bicgstab diverging ----------'
          write(*,*) '--------------- stopping ---------------'
          stop
        endif

        betatop = 0.d0
        do i = 1,nxny
          betatop = betatop + rhat(i) * r(i)
        enddo

        if (abs(alphatop) .lt. eps .or. abs(omega) .lt. eps ) then 
          beta_cg = 0.d0
          stop 'ach - bicgstab'
        else
          beta_cg = (betatop / alphatop) * (alpha_cg / omega)
        endif

        do i = 1,nxny
          p(i) = r(i) + beta_cg * ( p(i) - omega * w(i))
        enddo

        alphatop = betatop

      enddo

      write(*,*) 'residual: ', rmax
      write(*,*) 'bicgstab; max iterations reached; kmax = ',kmax
      stop

      ierr = 1

5     continue

      kfinal = k


c check that solution really satisfies linear system!!
      call dmat_vec(ws,u,
     &               ew_rc,ns_rc,
     &               nx,ny,dx,dy,dt,viscosity, neum,conv)

      rmax = 0.d0
      r2norm = 0.d0

      do i = 1,nxny
        r(i) = rhs(i) - ws(i)
        if ( abs(r(i)) .gt. rmax)  rmax = abs(r(i))
        r2norm = r2norm + r(i)**2.d0
      enddo
        r2norm = r2norm**0.5d0

      write(*,*) 'after solve, rmax = ',rmax
c      if (rmax .gt. tol*rhs_max) write(*,*) 'bicg tol not satisfied'
c      if (0.1d0*r2norm .gt. tol*rhs2norm) write(*,*) 'bicg tol not satf'
      if (0.1*rmax .gt. tol*rhs_max) stop 'bicg tol not satisfied'

      return
      end
