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

      subroutine update_navier_stokes(t, u, v, p, div,
     &            u_bc_ew_old, u_bc_ns_old, v_bc_ew_old, v_bc_ns_old,
     &            u_bc_ew,     u_bc_ns,     v_bc_ew,     v_bc_ns,
     &            u_bc_e_conv, u_bc_n_conv, v_bc_e_conv, v_bc_n_conv,
     &            dp_bc_ew,      dp_bc_ns, 
     &            b_old, b,
     &            An, A_ew, A_ns,
     &            tol_bicg, tol_cg,
     &            dx, dy, dt, viscosity, density, mcor,
     &            nx, ny, k, zero_dpdn, cond, neum, relative_tol, conv)

      implicit double precision(a-h,o-z)

      double precision us  (nx, ny)
      double precision vs  (nx, ny)
      double precision ps  (nx, ny)
      double precision ones(nx, ny)
c
c      double precision tu(nx, ny)
c      double precision tv(nx, ny)
c      double precision tp(nx, ny)
c      double precision error(nx,ny), error_p(nx,ny)  ! velc and pressure error

! solution at time t 
      double precision u   (nx, ny)
      double precision v   (nx, ny)
      double precision p   (nx, ny)
      double precision div (nx, ny)

      double precision u_no_pres(nx, ny)
      double precision v_no_pres(nx, ny)

! rhs terms
      double precision u_rhs   (nx, ny)
      double precision v_rhs   (nx, ny)
      double precision p_rhs   (nx, ny)

      double precision tol_cg(mcor)

! velocity boundary conditions
      double precision u_bc_ew_old(ny,2), u_bc_ns_old(nx,2)
      double precision v_bc_ew_old(ny,2), v_bc_ns_old(nx,2)

      double precision u_bc_ew(ny,2), u_bc_ns(nx,2)
      double precision v_bc_ew(ny,2), v_bc_ns(nx,2)

! vel bc's for div calc in neum case
      double precision u_bc_ew_neum(ny,2)
      double precision v_bc_ns_neum(nx,2)

! convective bcs; only used if conv .eq. .true.
      double precision u_bc_e_conv(ny),  u_bc_n_conv(nx)
      double precision v_bc_e_conv(ny),  v_bc_n_conv(nx)


! pressure boundary conditions
      double precision dp_bc_ew    (ny,2), dp_bc_ns    (nx,2) ! t = t_{n+1}


! Rhie-Chow interpolation terms (will also contain BC data)
! linear interpolation is used now
      double precision ew_rc (nx+1, ny)
      double precision ns_rc (nx , ny+1)

      double precision b_old (2*nx + 2*ny - 4, 2)
      double precision b(2*nx + 2*ny - 4, 2)

! PISO arrays
      double precision An(nx,ny) ! A^n times ones 
      double precision A_ew(nx-1, ny  ) ! interp (a_tilde)^-1 to e-w
      double precision A_ns(nx,   ny-1) ! interp (a_tilde)^-1 to n-s
      double precision A_e(ny) ! (a_tilde_p)^-1 along east
      double precision A_n(nx) ! (a_tilde_p)^-1 along north 
  
      double precision Anun(nx,ny), Anvn(nx,ny) ! A^n times vel^*
      double precision Hnun(nx,ny), Hnvn(nx,ny) ! H^n times vel^n
      double precision Hnus(nx,ny), Hnvs(nx,ny) ! H^n times vel^n
      
      double precision fs(nx,ny,2), f_old(nx,ny,2) ! forcing at vel^** ...

      double precision ew(nx-1,ny), ns(nx,ny-1)

      double precision dpdx(nx,ny), dpdy(nx,ny)

      double precision dpdx_edge(nx+1,ny), dpdy_edge(nx,ny+1)

      double precision one_v(nx,ny), col(nx,ny) ! for cg_solve matrix analysis
      

      character*4 filename

      logical rc_interp, zero_dpdn
      logical cond  ! preconditioned CG solver: yes/no
      logical neum ! Neumann bcs (North and East)
      logical conv ! convective bcs (North and East)
      logical relative_tol ! tol relative to rhs in cg_solve 

      logical cg_check ! true for recording mat,rhs,and x

      external vel_Hmat_vec, p_mat_vec, cg_solve_ns
      external filename, bicgstab
      external Hbar_vec, A_vec
      external cg_cond
      external u3_solve, l3_solve
      external get_p_diag
      external divergence_rc 
c      external calc_error

c for cg check:
      open(unit=24,file='matrix.dat',status='unknown')
      open(unit=25,file='rhs.dat',status='unknown')
      open(unit=26,file='x.dat',status='unknown')
 
      open(unit=44,file='ps.dat',status='unknown')
      open(unit=45,file='bc_ns.dat',status='unknown')
      open(unit=46,file='div_proj_after_intime.dat',status='unknown')
      open(unit=47,file='div_neum_after_intime.dat',status='unknown')
c end for cg check

      one = 1.d0

      integer_one = 1
      integer_two = 2

      call initialize(one, ones, nx, ny)

 
      call calc_b(b, 
     &            u_bc_ew, u_bc_ns,
     &            v_bc_ew, v_bc_ns,
     &            u_bc_e_conv, u_bc_n_conv,
     &            v_bc_e_conv, v_bc_n_conv,
     &            nx, ny, dx, dy, viscosity, neum, conv)


      call dp_edges(dpdx_edge, dpdy_edge, p, 
     &              dp_bc_ew, dp_bc_ns,
     &              nx, ny, dx, dy)


      call dp_centers(dpdx, dpdy, dpdx_edge, dpdy_edge,
     &                  nx, ny, dx, dy)

      rc_interp = .true.
! necessary data not available on first step
      if (k .eq. 1) rc_interp = .false.  

      call rhie_chow(u, v, p, ew_rc, ns_rc,
     &               dpdx, dpdy,
     &               dpdx_edge, dpdy_edge, 
     &               An, A_ew, A_ns,
     &               u_bc_ew_old, u_bc_ns_old, 
     &               v_bc_ew_old, v_bc_ns_old,
     &               dx, dy, dt, density, 
     &               nx, ny, rc_interp, neum) !  rc's at  n 
 
      call forcing(f_old, u, v, nx, ny, dx, dy, dt, t)

      call Hbar_vec(Hnun, u, ew_rc, ns_rc,
     &                nx, ny, dx, dy, dt, viscosity)

      call Hbar_vec(Hnvn, v, ew_rc, ns_rc,
     &                nx, ny, dx, dy, dt, viscosity) 


      call A_vec   (Anun, u, ew_rc, ns_rc,
     &              nx, ny, dx, dy, dt, viscosity, neum,conv) 
     
      call A_vec   (Anvn, v, ew_rc, ns_rc,
     &              nx, ny, dx, dy, dt, viscosity, neum,conv) 

      call A_vec   (An, ones, ew_rc, ns_rc,
     &              nx, ny, dx, dy, dt, viscosity, neum,conv) 


      call vel_rhs(u_rhs, u, Hnun, Anun, b_old, b,
     &             f_old, dpdx,
     &             nx, ny, dx, dy, dt, integer_one, 
     &             density)


      call vel_rhs(v_rhs, v, Hnvn, Anvn, b_old, b,
     &             f_old, dpdy,
     &             nx, ny, dx, dy, dt, integer_two, 
     &             density)


c-------initialize first quess for bicg--------------
      do i = 1,nx
        do j = 1,ny
          us(i,j) = u(i,j)
          vs(i,j) = v(i,j)
        enddo
      enddo

c-----end-initialize first guess for bicg------------     

      call bicgstab(us, u_rhs,
     &             ew_rc, ns_rc,
     &             viscosity, 
     &             dx, dy, dt,nx,ny,vel_Hmat_vec,ierr,kfinal_u,tol_bicg,
     &             neum,conv)


      call bicgstab(vs, v_rhs, 
     &             ew_rc, ns_rc,
     &             viscosity,   
     &             dx, dy, dt,nx,ny,vel_Hmat_vec,ierr,kfinal_v,tol_bicg,
     &             neum,conv)

c---------------------debug-------------------------------
c  begin check bicg
c
c      zero = 0.
c      do i=1,nx
c        do j=1,ny
c          write(25,*)  u_rhs(i,j)
c         write(26,*)  us(i,j)
c        enddo
c       enddo


c      do i =1,nx
c
c      do j =1,ny
c          call initialize(zero,one_v,nx,ny)
c          one_v(i,j)=1.d0
c        call vel_Hmat_vec(col,one_v,
c     &                    ew_rc,ns_rc,
c     &                    nx,ny,dx,dy,dt,viscosity,neum,conv)
c
c         do n =1,nx
c
c            do mm=1,ny
c             write(24,*) col(n,mm)
c            enddo
c         enddo
c
c        enddo
c       enddo
c       write(*,*) 'end bicg check'
c
c        stop
ccend check bicg
c----------------------end-debug----------------------------

      write(*,*) 'bicg iters (u,v):',kfinal_u, kfinal_v

c -------------------------debug----------------
c      
c          call matlabvtk(us,nx,ny,'us_n'//filename(k))
c          call matlabvtk(vs,nx,ny,'vs_n'//filename(k))
c          
c          stop
c-------------------------end-debug--------------------

 

c     call divergence(div, us, vs,
c    &                u_bc_ew, 
c    &                v_bc_ns,
c    &                div_rms, div_max,
c    &                dx, dy, nx, ny)

c     write(*,*) 'pred div (rms,max) ',div_rms,div_max

      write(*,*) '------------ end predictor -----------------'
      write(*,*) '  '

      do j = 1,ny 
        do i =1,nx
          ps(i,j) = p(i,j) 
        enddo
      enddo

      do m = 1, mcor

c-----------------------------------------------------------------------
c calculate pressure bc based on navier-stokes, if dpdx is not presc.
            call get_pbc_navier_stokes(
     &                   us, vs,
     &                   u_bc_ew_old, u_bc_ns_old, 
     &                   v_bc_ew_old, v_bc_ns_old,
     &                   u_bc_ew, u_bc_ns, v_bc_ew, v_bc_ns,
     &                   dp_bc_ew, dp_bc_ns, 
     &                   dt, nx, ny, 
     &                   dx, dy, density, viscosity, zero_dpdn)

c-----------------------------------------------------------------------

! Hbar^n(u^*)  and  Hbar^n(v^*)
          call Hbar_vec(Hnus, us, ew_rc, ns_rc,
     &                  nx, ny, dx, dy, dt, viscosity)
          call Hbar_vec(Hnvs, vs, ew_rc, ns_rc,
     &                  nx, ny, dx, dy, dt, viscosity) 

! forcing of u^* and v^*
          call forcing(fs, us, vs, nx, ny, dx, dy, dt, t) 

          call beta(u_no_pres, ew, ns, u, An, Hnus, Hnun, 
     &                b_old, b, fs, f_old,
     &                nx, ny, dx, dy, dt, density, integer_one)

          call beta(v_no_pres, ew, ns, v, An, Hnvs, Hnvn, 
     &                b_old, b, fs, f_old,
     &                nx, ny, dx, dy, dt, density, integer_two)

          call int_a_tilde(A_ew, A_ns, An, dx, dy, dt, nx, ny,
     &                     A_e, A_n)


          call pres_rhs(p_rhs, ew, ns, u_no_pres, v_no_pres,
     &                  u_bc_ew, u_bc_ns,
     &                  v_bc_ew, v_bc_ns,
     &                  u_bc_e_conv, u_bc_n_conv,
     &                  v_bc_e_conv, v_bc_n_conv,
     &                  nx, ny, dx, dy, dt, density, neum,conv)


c-----------------------debug-----------------------------------
c check p_rhs max
c          rhs_max = 0.d0
c          do i =1,nx
c            do j = 1,ny
c            if ( abs(p_rhs(i,j)) .gt. rhs_max) rhs_max = abs(p_rhs(i,j))
c            enddo
c           enddo
c           write(*,*) 'it',k,'cor',m,'rhs_max',rhs_max



c          if (k .eq. 1) then 
c          do i = 1,nx
c             write(*,*) 'p_rhs',p_rhs(nx,i), p_rhs(i,nx)
c             write(*,*)'no_pres', u_no_pres(nx,i), v_no_pres(i,nx)
c          enddo
c          stop
c          endif
c-----------------------end-debug----------------------------------

          if (neum) then

c  subtract ave from rhs 

      call sub_ave(p_rhs,nx,ny,.true.) !false to add

      write(*,*) tol_cg(1)
      write(*,*) tol_cg(m)

      call cg_solve_ns(ps, p_rhs, A_ew, A_ns, A_e, A_n,
     &                  dx, dy, dt, nx, ny, density, p_mat_vec, ierr,
     &                  kfinal, tol_cg(m), neum, relative_tol)

         else if (conv) then

c      call sub_ave(ps,nx,ny,.true.) !false to add

      call cg_solve_ns(ps, p_rhs, A_ew, A_ns, A_e, A_n,
     &                  dx, dy, dt, nx, ny, density, p_mat_vec, ierr,
     &                  kfinal, tol_cg(m), neum, relative_tol)

c---------------------------debug---------------------------------- 
c check p_mat_vec
c           if (k .eq. 3 .and. m .eq. 2) then
c
c            do i =1,nx
c              do j = 1,ny
c              p(i,j) = 1.d0
c              enddo
c            enddo
c            call p_mat_vec(p,p,A_ew,A_ns,A_e,A_n,
c     &                      nx,ny,dx,dy,dt,density,neum)
c
c            call matlabvtk(ps,nx,ny,'p_1')
c           write(*,*) 'check'
c          stop
c           endif
 
             
c            if (k .eq. 1) then
c             call p_mat_vec(ps,ps,A_ew,A_ns,A_e,A_n,
c     &                      nx,ny,dx,dy,dt,density,neum)
c            call matlabvtk(p_rhs,nx,ny,'p_1')
c            stop
c            endif
c-------------------------------end-debug-----------------------------

          else !~neum  & ~conv
     
            if (cond)  then

c          call sub_ave(p_rhs,nx,ny,.true.) !false to add
          call cg_cond(ps, p_rhs, A_ew, A_ns, A_e, A_n,
     &                 dx, dy, dt, nx, ny,
     &                 density, p_mat_vec,
     &                 ierr, kfinal, tol_cg(m), neum)

             else 


c          call sub_ave(p_rhs,nx,ny,.true.) !false to add
          call cg_solve_ns(ps, p_rhs, A_ew, A_ns, A_e, A_n,
     &                  dx, dy, dt, nx, ny, density, p_mat_vec, ierr, 
     &                  kfinal, tol_cg(m), neum, relative_tol)


             endif
 
          endif ! neum

c----------------------------------debug-------------------------- 
c begin check cg
      cg_check = .false.
      if (cg_check) then
c
      if (t .ge. 7 .and. m .eq. 2) then

      zero = 0.
      do i=1,nx
        do j=1,ny
          write(25,*)  p_rhs(i,j)
         write(26,*)  ps(i,j)
        enddo
       enddo


      do i =1,nx

      do j =1,ny
          call initialize(zero,one_v,nx,ny)
          one_v(i,j)=1.d0
        call p_mat_vec(col,one_v,A_ew,A_ns,A_e,A_n,
     &                      nx,ny,dx,dy,dt,density,.false.)

         do n =1,nx

            do mm=1,ny
             write(24,*) col(n,mm)
            enddo
         enddo

        enddo
       enddo
       write(*,*) 'kfinal',kfinal

        stop
        endif
      endif !cg_check
c
c end cg check
c--------------------------------end-debug---------------------------


          write(*,*) 'cg iter # for p:', 'corrector',m,kfinal

          call dp_edges(dpdx_edge, dpdy_edge, ps, 
     &                  dp_bc_ew, dp_bc_ns,
     &                  nx, ny, dx, dy)


          call dp_centers(dpdx, dpdy, dpdx_edge, dpdy_edge,
     &                    nx, ny, dx, dy)

          

          call cor_vel_un(us, u_no_pres, 
     &                    An, dpdx, 
     &                    density, nx, ny, dx, dy, dt)

          call cor_vel_un(vs, v_no_pres,  
     &                    An, dpdy, 
     &                    density, nx, ny, dx, dy, dt)


c-----debug---------------------
c
c        if (m .eq. 1) then
c          call matlabvtk(ps,nx,ny,'ps_n'//filename(k))
c          call matlabvtk(dpdx,nx,ny,'dpdx1_n'//filename(k))
c          call matlabvtk(dpdy,nx,ny,'dpdy1_n'//filename(k))
c          call matlabvtk(dpdx_edge,nx+1,ny,'dpdx1_edge_n'//filename(k))
c          call matlabvtk(dpdy_edge,nx,ny+1,'dpdy1_edge_n'//filename(k))
c
c          elseif (m .eq. 2) then
c          call matlabvtk(ps,nx,ny,'p_n'//filename(k))
c          call matlabvtk(dpdx,nx,ny,'dpdx2_n'//filename(k))
c          call matlabvtk(dpdy,nx,ny,'dpdy2_n'//filename(k))
c          call matlabvtk(dpdx_edge,nx+1,ny,'dpdx2_edge_n'//filename(k))
c          call matlabvtk(dpdy_edge,nx,ny+1,'dpdy2_edge_n'//filename(k))
c        endif


c-----end-debug-------------------


c-------------------------------------debug----------------------
c         call divergence_rc(div, us, vs,
c     &                    u_bc_ew, v_bc_ns,
c     &                    dpdx, dpdy,
c     &                    dpdx_edge, dpdy_edge, 
c     &                    An, A_ew, A_ns,
c     &                    div_rms, div_max, density,
c     &                    dx, dy, dt, nx, ny)

c         call divergence(div, us, vs,
c     &                    u_bc_ew,
c     &                    v_bc_ns,
c     &                    div_rms, div_max,
c     &                    dx,dy,nx,ny)
 

c       call matlabvtk(div,nx,ny,'dive')
c      stop


c         write(*,*) 
c         write(*,*) 'corrected div (rms,max)',div_rms, div_max
c         write(*,*) 
c         ivtk_error_track = ivtk_error_track + 1

c         write(20,*)  

c         write(71,*) m, div_max
c--------------------------------------end-debug-----------------------
   
          write(*,*) ' -------------- end of corrector ',m,'----------'


          write(135,*) u((nx-1)/2,(ny-1)/2),
     &                 us((nx-1)/2,(ny-1)/2)


        enddo  ! end of correction loop


c       write(198,*) t, div_rms
c       write(199,*) t, div_max 
c       call flush(198)
c       call flush(199)


c--------find  p^(n+1), u^(n+1), v^(n+1) --------------------
        do j = 1,ny
          do i= 1,nx
            p(i,j) = ps(i,j) 
            u(i,j) = us(i,j)
            v(i,j) = vs(i,j)
          enddo
        enddo

        do i = 1, 2*nx + 2*ny - 4
          b_old(i,1) = b(i,1)
          b_old(i,2) = b(i,2)
        enddo

        do j = 1, 2
          do i = 1, ny
            u_bc_ew_old(i,j) = u_bc_ew(i,j)
            v_bc_ew_old(i,j) = v_bc_ew(i,j) 
          enddo
        enddo
        do j = 1, 2
          do i = 1, nx
            u_bc_ns_old(i,j)  = u_bc_ns(i,j)
            v_bc_ns_old(i,j)  = v_bc_ns(i,j) 
          enddo
        enddo

c----debug: check divergence in neum case----------------
c          call divergence_rc(div,u,v,
c     &                       u_bc_ew, v_bc_ns,
c     &                       dpdx, dpdy,
c     &                       dpdx_edge, dpdy_edge,
c     &                       An, A_ew, A_ns,
c     &                       div_rms, div_max, density,
c     &                       dx, dy, dt, nx, ny)
c 
c      write(45,*) t,div_rms, div_max
 
c----end debug: check divergence in neum case----------------

c------debug: calc divergence after NS step---------------
       if (.not. neum) then
         call divergence_rc(div, u, v,
     &                    u_bc_ew, v_bc_ns,
     &                    dpdx, dpdy,
     &                    dpdx_edge, dpdy_edge, 
     &                    An, A_ew, A_ns,
     &                    div_rms, div_max, density,
     &                    dx, dy, dt, nx, ny)

c         call divergence(div, u, v,
c     &                    u_bc_ew,
c     &                    v_bc_ns,
c     &                    div_rms, div_max,
c     &                    dx,dy,nx,ny)
c
c        call matlabvtk(div,nx,ny,'div_proj_after'//filename(k))
         
         write(46,*) k,div_max,div_rms

        endif !proj

        if (neum) then

          do j = 1,ny !for neum div calc
             u_bc_ew_neum(j,1) = u_bc_ew(j,1)
             u_bc_ew_neum(j,2) = u(nx,j)
          enddo

          do i = 1,nx ! for neum div calc 
             v_bc_ns_neum(i,1) = v_bc_ns(i,1)
             v_bc_ns_neum(i,2) = v(i,ny)
          enddo
       
          call divergence_rc(div, u, v,
     &                    u_bc_ew_neum, v_bc_ns_neum,
     &                    dpdx, dpdy,
     &                    dpdx_edge, dpdy_edge, 
     &                    An, A_ew, A_ns,
     &                    div_rms, div_max, density,
     &                    dx, dy, dt, nx, ny)

c          call divergence(div, u, v,
c     &                    u_bc_ew_neum,
c     &                    v_bc_ns_neum,
c     &                    div_rms, div_max,
c     &                    dx,dy,nx,ny)

c        call matlabvtk(div,nx,ny,'div_neum_after'//filename(k))

         write(47,*) k,div_max,div_rms
 
 
       endif !neum

c------debug: calc divergence after NS step---------------

        return
        end      


