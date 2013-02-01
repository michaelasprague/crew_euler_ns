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

      subroutine cg_solve_ns(phi, rhs, A_ew, A_ns, A_e, A_n,
     &                    dx, dy, dt, nx, ny,
     &                    density, dmat_vec, ierr, kfinal, tol, neum,
     &                    relative_tol)

      implicit double precision (a-h,o-z)

C primary variables
c      double precision u(nx*ny), v(nx*ny), p(nx*ny)
      double precision phi(nx*ny)
      double precision A_ew((nx-1)*ny)
      double precision A_ns(nx*(ny-1))
      double precision A_e(ny)
      double precision A_n(nx)

C RHS forcing in pressure equation
      double precision rhs(nx*ny)

C   Conjugate gradient arrays
      double precision w(nx*ny), r(nx*ny), p(nx*ny)

      logical neum
      logical fix_mid 
      logical relative_tol

      external dmat_vec

      tol_cg = tol

      nxny = nx*ny

      kmax = nxny*100

      kfinal = 0
      k = 0

      ierr = 0

      eps = 1.d-40

c mas
      istart = 2 ! fixing first point (istart=2)
      phi(1)  = 0.d0
      r(1)  = 0.d0
      p(1)  = 0.d0

c      if (neum) istart = 1


      
      iend = nxny ! fixing last point (iend=nxny-1)
c      phi(nxny) = 0.d0
c      r(nxny) = 0.d0
c      p(nxny) = 0.d0

c fixing middle point (fix_mid=true, istart=1, iend =nxny)

      fix_mid = .false.
      mid = nxny/2
c      phi(mid) = 0.d0
c      r(mid) = 0.d0
c     p(mid) = 0.d0
 

c mas
c     w = A * phi
      call dmat_vec(w, phi, A_ew, A_ns, A_e, A_n,
     &               nx, ny, dx, dy, dt, density, neum)

      rmax = 0.d0
      rhs_max = 0.d0
      r2norm = 0.d0
      rhs2norm = 0.d0

      do i = istart,iend
        
        if (fix_mid .and. i .eq. mid) then
          r(i) = 0.d0
          p(i) = 0.d0
        else

          r(i) = rhs(i) - w(i)
          p(i) = r(i)
          if ( abs(r(i)) .gt. rmax)  rmax = abs(r(i))
          if ( abs(rhs(i)) .gt. rhs_max)  rhs_max = abs(rhs(i))
          r2norm = r2norm + r(i)**2.d0
          rhs2norm = rhs2norm + rhs(i)**2.d0

        endif
      enddo
        
        if (.not. relative_tol)  then
          rhs_max = 1.d0
        endif

        tol_cg = tol_cg * rhs_max

        if (abs(tol_cg) .lt. 1.d-15) then 
         tol_cg = 1.d-15
        endif
       
c       write(*,*) 'IN CG:','rhs_max', rhs_max, 'tol,tol_cg',tol,tol_cg

        r2norm = r2norm**0.5d0
        rhs2norm = rhs2norm**0.5d0


      if (rmax .lt. tol_cg)  goto 5

      alphatop = 0.d0
      do i = 1,nxny
        alphatop = alphatop + r(i) * r(i)
      enddo

      do k = 1,kmax
         
c  mas
c       w = A * p
        call dmat_vec(w,p,A_ew,A_ns, A_e, A_n,
     &               nx,ny,dx,dy,dt, density, neum)

        alphabot = 0.d0
        do i = 1,nxny
          alphabot = alphabot + p(i) * w(i)
        enddo

        if (abs(alphabot) .lt. eps .and. alphatop .gt. eps) then
          write(*,*) 'cg_solve_ns:alphabot  is', alphabot 
          write(*,*) 'cg_solve_ns:alphatop  is', alphatop
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
        r2norm = 0.d0

        do i = istart,iend
          if (fix_mid .and. i .eq. mid) then
           r(i) = 0.d0
           phi(i) = 0.d0
          else
 
           phi(i) = phi(i) + alpha_cg * p(i)

           r(i) = r(i) - alpha_cg*w(i)  ! r^(k+1) = r^(k) - alpha A p^(k)

          if ( abs(r(i)) .gt. rmax)  rmax = abs(r(i))
          r2norm = r2norm + r(i)**2.d0
         endif

        enddo
          r2norm = r2norm**0.5d0
 
c        write(*,*) phi(4*41+5)

c uncomment if you want to try to fix the pressure
c       phi(1) = 0.d0
        
        if (rmax .lt. tol_cg)  goto 5

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

        do i = istart,iend
          if (fix_mid .and. i .eq. mid) then
           p(i) = 0.d0
          else

          p(i) = r(i) + beta_cg * p(i)
          endif

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

      kfinal = k

5     continue


c     ave = 0.d0
c     do i = 1, nxny
c       ave = ave + phi(i)
c     enddo
c     ave = ave / float(nxny)
c     do i = 1, nxny
c       phi(i) = phi(i) - ave
c     enddo

c check that solution really satisfies linear system!!
      call dmat_vec(w, phi, A_ew, A_ns, A_e, A_n,
     &             nx, ny, dx, dy, dt, density, neum)

      rmax = 0.d0
      r2norm = 0.d0
      do i = istart,iend
       if (fix_mid .and. i .eq. mid) then
         r(i) = 0.d0
       else

        r(i) = rhs(i) - w(i)
        if ( abs(r(i)) .gt. rmax)  rmax = abs(r(i))
        r2norm = r2norm + r(i)**2.d0
       endif

      enddo
      r2norm = sqrt(r2norm)

      write(*,*) 'after solve, rmax = ',rmax
      if (0.1*rmax .gt. tol_cg) stop 'cg tolerance not satisfied'

      kfinal = k

      return
      end
