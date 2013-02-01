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

      subroutine proj_cg_solve(phi, rhs, 
     &                    dx, dy, nx, ny,
     &                    dmat_vec, ierr, kfinal, tol)

      implicit double precision (a-h,o-z)

C primary variables
c      double precision u(nx*ny), v(nx*ny), p(nx*ny)
      double precision phi(nx*ny)

      double precision u(nx*ny), v(nx*ny)
      
      double precision bc_ew(ny,2), bc_ns(nx,2)

C RHS forcing in pressure equation
      double precision rhs(nx*ny)

C   Conjugate gradient arrays
      double precision w(nx*ny), r(nx*ny), p(nx*ny)

      external dmat_vec

      nxny = nx*ny

      kmax = nxny*1000

      kfinal = 0
      k = 0

      ierr = 0

      eps = 1.d-30

c mas
      phi(1)  = 0.d0
      r(1)  = 0.d0
      p(1)  = 0.d0
      istart = 1

      call dmat_vec(w, phi, nx, ny, dx, dy)

      rmax = 0.d0
      do i = 1,nxny
        r(i) = rhs(i) - w(i)
        p(i) = r(i)
        if ( abs(r(i)) .gt. rmax)  rmax = abs(r(i))
      enddo
      if (rmax .lt. tol)  goto 5

      alphatop = 0.d0
      do i = istart,nxny
        alphatop = alphatop + r(i) * r(i)
      enddo

      do k = 1,kmax

c  mas
c       p(1) = 0.d0
        call dmat_vec(w,p,nx,ny,dx,dy)

        alphabot = 0.d0
        do i = istart,nxny
          alphabot = alphabot + p(i) * w(i)
        enddo

        if (abs(alphabot) .lt. eps .and. alphatop .gt. eps) then
          write(*,*) 'alphabot  is', alphabot 
          write(*,*) 'alphatop  is', alphatop
          write(*,*) 'stopping'
          stop
        endif

        if (abs(alphabot) .lt. eps) then
          alpha_cg = 0.d0
          stop 'check this -- cgsolve'
        else
          alpha_cg = alphatop / alphabot
        endif

        rmax = 0.d0
        do i = istart, nxny

          phi(i) = phi(i) + alpha_cg * p(i)

          r(i) = r(i) - alpha_cg*w(i)  ! r^(k+1) = r^(k) - alpha A p^(k)

          if ( abs(r(i)) .gt. rmax)  rmax = abs(r(i))

        enddo
     
c        write(*,*) phi(4*41+5)

c uncomment if you want to try to fix the pressure
c       phi(1) = 0.d0

        if (rmax .lt. tol)  goto 5

        if (rmax .gt. 1d6) then
          write(*,*) 'k = ',k
          write(*,*) '------------- cg diverging -------------'
          write(*,*) '--------------- stopping ---------------'
          stop
        endif

        betatop = 0.d0
        do i = 1,nxny
          betatop = betatop + r(i) * r(i)
        enddo
      
        beta_cg = betatop / alphatop

        do i = istart,nxny
          p(i) = r(i) + beta_cg * p(i)
        enddo

        alphatop = betatop

      enddo

      write(*,*) 'max iterations reached; kmax = ',kmax
      write(*,*) 'solution written to fort.71'
      do i = 1, nxny
        write(71,*) i, p(i)
      enddo

      stop
      ierr = 1

5     continue

      kfinal = k

c     ave = 0.d0
c     do i = 1, nxny
c       ave = ave + phi(i)
c     enddo
c     ave = ave / float(nxny)
c     do i = 1, nxny
c       phi(i) = phi(i) - ave
c     enddo


c check that solution really satisfies linear system!!
      call dmat_vec(w,phi,nx,ny,dx,dy)

      rmax = 0.d0
      do i = istart,nxny
        r(i) = rhs(i) - w(i)
        if ( abs(r(i)) .gt. rmax)  rmax = abs(r(i))
      enddo
!     write(*,*) 'after solve, rmax = ',rmax
      if (0.1*rmax .gt. tol) stop 'cg tolerance not satisfied'


      return
      end
