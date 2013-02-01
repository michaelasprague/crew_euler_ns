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

      program driver

c solving for pressure directy

c 2-D incompressible Navier-Stokes solver
c employs cell-centered finite-volumes with Rhie-Chow interpolation
c CG Solver for pressure
c PISO corector predictor w/ trap. time stepping
c
c restricted to a rectangular grid with rectangular elements
c
c Coding started on 12 January 2011
c
c contact: michael.a.sprague@nrel.gov, michaelasprague@gmail.com
c

      implicit double precision (a-h,o-z)

      parameter (nx = 31)
      parameter (ny = 31) 
      parameter (mcor = 2)

      double precision us  (nx, ny)
      double precision vs  (nx, ny)
      double precision ps  (nx, ny)
      double precision p1  (nx, ny)

! solution at time t 
      double precision u  (nx, ny)
      double precision v  (nx, ny)
      double precision p  (nx, ny)
      double precision pold(nx, ny)
      double precision div(nx, ny)

! rhs terms
      double precision u_rhs   (nx, ny)
      double precision v_rhs   (nx, ny)
      double precision p_rhs   (nx, ny)

! velocity boundary conditions
      double precision u_bc_ew_n(ny,2), u_bc_ns_n(nx,2)  ! t = t_n
      double precision v_bc_ew_n(ny,2), v_bc_ns_n(nx,2)
      double precision u_bc_ew_n1(ny,2), u_bc_ns_n1(nx,2) ! t = t_{n+1}
      double precision v_bc_ew_n1(ny,2), v_bc_ns_n1(nx,2)

! convective bcs
      double precision u_bc_e_conv(ny),  u_bc_n_conv(nx)
      double precision v_bc_e_conv(ny),  v_bc_n_conv(nx)

      double precision dp_bc_ew_n(ny,2), dp_bc_ns_n(nx,2) ! t = t_{n+1}
      double precision dp_bc_ew_n1(ny,2), dp_bc_ns_n1(nx,2) ! t = t_{n+1}
      double precision dpdx(nx,ny), dpdy(nx,ny)
      double precision dpdx_edge(nx+1,ny), dpdy_edge(nx,ny+1)

! b terms in momentum predictor eq.
      ! second coordinate is for velocity: 1 for u; 2 for v

      double precision bn( 2*nx + 2*ny - 4, 2)
      double precision bn1(2*nx + 2*ny - 4, 2)

!  Taylor Green arrays
      double precision tu(nx,ny), tv(nx,ny), tp(nx,ny)
      double precision error(nx,ny), error_p(nx,ny)  ! velc and pressure error

! Rhie-Chow interpolation terms (will also contain BC data)
! linear interpolation is used now

      double precision ew_rc    (nx+1, ny)
      double precision ns_rc    (nx, ny+1) 

      double precision tol_cg   (mcor)

! PISO arrays
      double precision An(nx,ny)     ! A^n times ones 
      double precision A_ew(nx-1,ny) ! interp (a_tilde)^-1 to e-w
      double precision A_ns(nx,ny-1) ! interp (a_tilde)^-1 to n-s
  
      double precision Anun(nx,ny), Anvn(nx,ny) ! A^n times vel^*
      double precision Hnun(nx,ny), Hnvn(nx,ny) ! H^n times vel^n
      double precision Hnus(nx,ny), Hnvs(nx,ny) ! H^n times vel^n
      
      double precision fs(nx,ny,2), fn(nx,ny,2) ! forcing at vel^** ...

      double precision ew(nx-1,ny), ns(nx,ny-1)
      
c vtk 

      integer ns_vtk_option ! calc mean or subtract given mean

      character*4 filename

      logical zero_dpdn
      logical cond  !preconditioned CG: on/off
      logical conv  ! convective bcs                                           
      logical neum  ! Neumann bcs (North and South)
      logical relative_tol ! tol relative to rhs in cg_solve 

      logical if_matlabvtk

      logical isnan

      logical if_kovas

      logical restart

      external vel_Hmat_vec, p_mat_vec, cg_solve_ns
      external filename,initialize, bicgstab
      external Hbar_vec, A_vec
      external cg_cond
      external u3_solve, l3_solve
      external get_p_diag
      external calc_error

      open(unit=19,file='u_center.dat',status='unknown')

      open(unit=17,file='restart_old.dat',status='unknown')
      open(unit=18,file='restart_new.dat',status='unknown')

      open(unit=19,file='u_center.dat',status='unknown')
c     open(unit=21,file='u_iters.dat',status='unknown')
c     open(unit=22,file='v_iters.dat',status='unknown')

c for cg check:
      open(unit=24,file='matrix.dat',status='unknown')
      open(unit=25,file='rhs.dat',status='unknown')
      open(unit=26,file='x.dat',status='unknown')
c end for cg check

      open(unit=30,file='rms.dat',status='unknown')
      open(unit=31,file='div.dat',status='unknown')

      open(unit=40,file='params.dat',status='unknown')

      open(unit=47,file='divergence.dat',status='unknown')

! ---------------------------------------------------------------------
! Problem parameters
! ---------------------------------------------------------------------
c c1

C Choose test case
      if_kovas = .true.

      if (if_kovas) then
        ! these parameters are set for the kovas case; don't modify
        if_matlabvtk = .false.
        zero_dpdn    = .true.
        cond         = .false. !preconditioned CG: on/off
        neum         = .false. !use Neumann bcs (East and North)
        relative_tol = .true. !relative tol to p_rhs in cg_solve
        restart      = .false. 

        dt_vtk_out = 0.1d0

        xmin = -0.5d0 
        xmax = +0.5d0 
        ymin = -0.5d0 
        ymax = +0.5d0 

        tmax = 3.d0
        cn = 0.5d0

        tol_cg(1)   = 1.d-10
        tol_cg(2)   = 1.d-10
        tol_bicg = 1.d-10

        viscosity = 0.025d0  ! Re = 40
        density   = 1.d0

        velmax = 2.62  ! from continuous form, 2.6191 
      else
        if_matlabvtk = .false.
        zero_dpdn    = .true.
        cond         = .false. !preconditioned CG: on/off
        neum         = .false. !use Neumann bcs (East and North)
        relative_tol = .true.  !relative tol to p_rhs in cg_solve
        restart      = .false. 

        dt_vtk_out = 0.1d0

        xmin = -0.5d0 
        xmax = +0.5d0 

        ymin = -0.5d0 
        ymax = +0.5d0 

        tmax = 1.d0
        cn = 0.5d0

        tol_cg(1)   = 1.d-10
        tol_cg(2)   = 1.d-10

        tol_bicg = 1.d-10

        tol_bicg_p = 1.d-8

        viscosity = 0.025d0  ! Re = 40
        density   = 1.d0

        velmax = 2.62  ! from continuous form, 2.6191 

      endif



      dx = (xmax - xmin) / float(nx)
      dy = (ymax - ymin) / float(ny)

      pi = acos(-1.d0)
      

      dx_min = sqrt(dx**2 + dy**2)
 
      dt = cn * dx_min / velmax

      i_t_max = int(tmax / dt)

      write(40,*) 'dt,dx,dy,cn ',dt,dx,dy,cn


      write(40,*) 'Peclet',(velmax*dx) / viscosity
 
! ---------------------------------------------------------------------
c c3  vtk control
 
      dt_vtk     = 0.d0
      i_vtk_out  = 0

      ns_vtk_option = 1 ! subtract given umean and vmean

c mean velocity for outputvtk_ns (subtracts it)
c also, initialize routine uses it

      umean = 0.d0
      vmean = 0.d0
      pmean = 0.d0
     

c this is how it looks down there in vtk  output
c       if (dt_vtk .ge. dt_vtk_out) then
c          dt_vtk = 0.d0
c          i_vtk_out = i_vtk_out + 1
c          call outputvtk(....i_vtk_out....
c        else
c          dt_vtk = dt_vtk + dt
c        endif


      t = 0.

! set everything to zero
      zero = 0.
      one = 1.d0
      almost_one = 1.1d0

      call initialize(umean,u,nx,ny)
      call initialize(vmean,v,nx,ny)
      call initialize(zero,p,nx,ny)
      call initialize(zero,pold,nx,ny)
      call initialize(zero,ps,nx,ny)

      call initialize(zero,error,nx,ny)

c initialize for Kovasznay 
       call kovas(u, v,  p,xmin,ymin,dx,dy,
     &            nx,ny,t,viscosity,density)
 
       call kovas(tu,tv,tp,xmin,ymin,dx,dy,
     &            nx,ny,t,viscosity,density)

      call calc_error(error,error_p,u,v,p,tu,tv,tp,
     &                u_err_max, u_err_min,
     &                rmsu,rmsp,nx,ny)
      write(*,*) 'initialization error:  rmsu, rmsp ',rmsu,rmsp


      small = 0.0
      do j = 1,ny
        do i = 1,nx
           scale = 1. + small*rand()
           u(i,j) = scale * u(i,j)
           scale = 1. + small*rand()
           v(i,j) = scale * v(i,j)
           scale = 1. + small*rand()
           p(i,j) = scale * p(i,j)

           us(i,j) = u(i,j)
           vs(i,j) = v(i,j)

           pold(i,j) = p(i,j)

        enddo
      enddo

c---------find b_0---------------

c inflow bc's (for choking the code)
c       flux_tot = 1.d-1
c       call inflow_bcs(u_bc_ew_n, u_bc_ns_n, v_bc_ew_n, v_bc_ns_n,
c     &                        dp_bc_ew, dp_bc_ns,
c     &                         nx, ny, dx, dy, viscosity,flux_tot)


c uniform bcs: uv=1 all around

c       call uniform_bcs(u_bc_ew_n, u_bc_ns_n, 
c    &                   v_bc_ew_n, v_bc_ns_n,
c    &                   nx, ny, dx, dy)
c 

c--------debug---------------------------
c       do i=1,nx
c         write(*,*) v_bc_ns_n(i,1),v_bc_ns_n(i,2) 
c       enddo
c       stop
c---------end-debug--------------------------     
           

      call get_kovas_pbc(u_bc_ew_n, u_bc_ns_n, 
     &                   v_bc_ew_n, v_bc_ns_n,
     &                   dp_bc_ew_n, dp_bc_ns_n,
     &                   xmin, xmax, ymin, ymax,
     &                   t, nx, ny, dx, dy, viscosity)

      call get_kovas_pbc(u_bc_ew_n1, u_bc_ns_n1, 
     &                   v_bc_ew_n1, v_bc_ns_n1,
     &                   dp_bc_ew_n1, dp_bc_ns_n1,
     &                   xmin, xmax, ymin, ymax,
     &                   t, nx, ny, dx, dy, viscosity)

      call calc_b(bn, 
     &          u_bc_ew_n, u_bc_ns_n,
     &          v_bc_ew_n, v_bc_ns_n,
     &          u_bc_ew_n, u_bc_ns_n,
     &          v_bc_ew_n, v_bc_ns_n,
     &          nx, ny, dx, dy, viscosity, neum, conv)

c calculate pressure bc based on navier-stokes, if dpdx is not presc.
c "true" for dpdn=0

        call get_pbc_navier_stokes(
     &               u, v, 
     &               u_bc_ew_n, u_bc_ns_n, v_bc_ew_n, v_bc_ns_n,
     &               u_bc_ew_n, u_bc_ns_n, v_bc_ew_n, v_bc_ns_n,
     &               dp_bc_ew_n1, dp_bc_ns_n1, 
     &               dt, nx, ny, 
     &               dx, dy, density, viscosity, zero_dpdn)

c---------start check divergence of BC's-------------
   
      s = 0.d0
      do  i=1,nx
        s = s + v_bc_ns_n(i,2)*dx - v_bc_ns_n(i,1)*dx
      enddo
      do  j=1,ny
        s = s + u_bc_ew_n(j,2)*dy - u_bc_ew_n(j,1)*dy
      enddo
      write(*,*) 's=',s

      if (abs(s) .gt. 1d-10) stop 'bcs not div free'

c---------end check divergence of BC's-------------

c-----------------------------------------------------------------------
c     start time integration
c-----------------------------------------------------------------------


      if (restart) then
        read(17,*) t, nxtmp, nytmp
        if (nxtmp .ne. nx .or. nytmp .ne. ny) stop 'restart not matched'
        do j = 1, ny
          do i = 1, nx
             read(17,*) u(i,j), v(i,j), p(i,j)
          enddo
        enddo
      endif

      write(19,*) t, u((nx+1)/2,(ny+1)/2)
      call flush(19)

c     call divergence_rc(div, u, v,
c    &                   u_bc_ew_n, v_bc_ns_n,
c    &                   dpdx, dpdy,
c    &                   dpdx_edge, dpdy_edge,
c    &                   An, A_ew, A_ns,
c    &                   div_rms, div_max,
c    &                   dx,dy,nx,ny)

c      call matlabvtk(u,nx,ny,'u_vel'//filename(ivtk_track))


      call calc_error(error,error_p,u,v,p,tu,tv,tp,
     &                u_err_max, u_err_min,
     &                rmsu,rmsp,nx,ny)

      call outputvtk_ns(u, v, p, div,
     &                  xmin, ymin, dx, dy,
     &                  nx,ny,'output'//filename(0),
     &                  umean,vmean,pmean,ns_vtk_option)

      call  outputvtk_bound(u_bc_ew_n, u_bc_ns_n,
     &                      v_bc_ew_n, v_bc_ns_n,
     &                      xmin, ymin, xmax, ymax,
     &                      umean, vmean,
     &                      dx,dy,nx,ny,
     &                      'output_bound'//filename(0))


c----------------------------------------------------------------------
c----------------- start testing zone ---------------------------------


c----------------- end testing zone ---------------------------------
c----------------------------------------------------------------------


      do i_t = 1, i_t_max  

        write(*,*) '-------------------------------------------'
        write(*,*) 'iteration ', i_t
        write(*,*) 100.*float(i_t) / float(i_t_max), '% complete'


C this needs to be modified if boundary conditions are time dependent 
c  -- this is really just here to get the pressure bcs
c       call get_kovas_pbc(u_bc_ew_n1, u_bc_ns_n1, 
c    &                     v_bc_ew_n1, v_bc_ns_n1,
c    &                     dp_bc_ew_n1, dp_bc_ns_n1,
c    &                     xmin, xmax, ymin, ymax,
c    &                     t+dt, nx, ny, dx, dy, viscosity)

c       call inflow_bcs(u_bc_ew_n1, u_bc_ns_n1, v_bc_ew_n1, v_bc_ns_n1,
c     &                        dp_bc_ew, dp_bc_ns,
c     &                         nx, ny, dx, dy, viscosity,flux_tot)
       
c       call uniform_bcs(u_bc_ew_n1, u_bc_ns_n1,
c    &                  v_bc_ew_n1, v_bc_ns_n1,
c    &                         nx, ny, dx, dy)

     
        call update_navier_stokes(t, u, v, p, div,
     &            u_bc_ew_n,  u_bc_ns_n,  v_bc_ew_n,  v_bc_ns_n,
     &            u_bc_ew_n1, u_bc_ns_n1, v_bc_ew_n1, v_bc_ns_n1,
     &            u_bc_e_conv, u_bc_n_conv, v_bc_e_conv, v_bc_n_conv,
     &            dp_bc_ew_n1, dp_bc_ns_n1, 
     &            bn, bn1,
     &            An, A_ew, A_ns,
     &            tol_bicg, tol_cg,
     &            dx, dy, dt, viscosity, density, mcor,
     &            nx, ny, i_t, zero_dpdn, cond, neum,
     &            relative_tol, conv)
           
c---------------------------------debug----------------------
c checking b terms

c      if (i_t .eq. 1) then
c        do i = 1, 2*nx+2*ny-4
c           write(*,*) bn(i,1),bn(i,2)
c        enddo
c      stop
c      endif
      
c check dp_bc

c      do i = 1,nx
c        write(*,*) dp_bc_ew_n1(i,2),dp_bc_ns_n1(i,2)
c      enddo
    
c--------------------------------end-debug---------------------


c--------Kovasznay flow-------------

c         call kovas(tu,tv,tp,xmin,ymin,dx,dy,
c     &             nx,ny,t+dt,viscosity,density)


c--------end Kovasznay flow-----------------


c---- Euclidean norm of the error-----------
        call calc_error(error, error_p, u, v, p, tu, tv, tp,
     &                  u_err_max, u_err_min,
     &                  rmsu, rmsp, nx, ny)

        write(*,*) 'u_err_max, u_err_min', u_err_max, u_err_min
        write(*,*) 'err rms (u^n+1,v^n+1), rms (p^n+1)  ',rmsu,rmsp

        write(30,*) t,  rmsu,    rmsp
        call flush(30)

c-----end Euclidean norm------------------


c this is how it looks up there
c      dt_vtk_out =0.1d0  or =0.0 to output every time step
c      dt_vtk = 0.d0
c      i_vtk_out = 0

       if (dt_vtk .ge. dt_vtk_out) then
c----------- START VTK ---------------
         write(*,*) 'START VTK'
          dt_vtk = 0.d0
          i_vtk_out = i_vtk_out + 1


          call dp_edges(dpdx_edge, dpdy_edge, p,
     &              dp_bc_ew_n1, dp_bc_ns_n1,
     &              nx, ny, dx, dy)


           call dp_centers(dpdx, dpdy, dpdx_edge, dpdy_edge,
     &                     nx, ny, dx, dy)


c          call divergence(div, u, v,
c     &                      u_bc_ew,
c     &                      v_bc_ns,
c     &                      div_rms, div_max,
c     &                      dx, dy, nx, ny)

         if (i_t .gt. 1) then

           call divergence_rc(div,u,v,
     &                       u_bc_ew_n1, v_bc_ns_n1,
     &                       dpdx, dpdy,
     &                       dpdx_edge, dpdy_edge,
     &                       An, A_ew, A_ns,
     &                       div_rms, div_max, density,
     &                       dx, dy, dt, nx, ny)
          

          endif

          write(*,*) 'it=', i_t
          write(*,*) 'div_rms', div_rms

          call calc_error(error,error_p,u,v,p,tu,tv,tp,
     &                    u_err_max, u_err_min,
     &                    rmsu,rmsp,nx,ny)

          call outputvtk_ns(u, v, p, div, xmin,ymin,dx,dy,nx,ny,
     &                       'output'//filename(i_vtk_out),
     &                       umean, vmean, pmean, ns_vtk_option)

          call  outputvtk_bound(u_bc_ew_n1, u_bc_ns_n1,
     &                          v_bc_ew_n1, v_bc_ns_n1,
     &                          xmin, ymin, xmax, ymax,
     &                          umean, vmean,
     &                          dx,dy,nx,ny,
     &                          'output_bound'//filename(i_vtk_out))

          write(*,*) 'END VTK'

c-------------matlab vtk----------------------
          if (if_matlabvtk) then
            call matlabvtk(u,nx,ny,'u_vel_n'//filename(i_vtk_out))
            call matlabvtk(v,nx,ny,'v_vel_n'//filename(i_vtk_out))
          endif
c          call matlabvtk(p,nx,ny,'p_n'//filename(i_vtk_out))


c-----------end-matlab vtk----------------------


          write(*,*)  '-------------------------------'
          !if (isnan(div_rms)) stop '"div_rms" is a NaN'
          !if (isnan(div_max)) stop '"div_max" is a NaN'
          write(*,*)  'div (rms,max):', div_rms, div_max
          write(*,*)  '-------------------------------'


        else
          dt_vtk = dt_vtk + dt
        endif

c-------------END VTK-------------------

        t = t + dt

        write(19,*) t, u((nx+1)/2,(ny+1)/2)
        call flush(19)
     
      enddo ! timeloop
c-----------------------------------------------------------------------
c     end   time integration
c-----------------------------------------------------------------------

c write restart file
      write(18,*) t, nx, ny
      do j = 1, ny
        do i = 1, nx
           write(18,*) u(i,j), v(i,j), p(i,j)
        enddo
      enddo

      call divergence_rc(div,u,v,
     &                   u_bc_ew_n1, v_bc_ns_n1,
     &                   dpdx, dpdy,
     &                   dpdx_edge, dpdy_edge,
     &                   An, A_ew, A_ns,
     &                   div_rms, div_max, density,
     &                   dx, dy, dt, nx, ny)
         
       write(*,*) 'div_rms, div_max:', div_rms, div_max

      write(*,*) '------------------------------------------------'
      write(*,*) 'u(1,1) :', u(1,1), v(1,1)
      write(*,*) 'u(2,1) :', u(2,1), v(2,1)
      write(*,*) 'u(1,2) :', u(1,2), v(1,2)
      write(*,*) '------------------------------------------------'
      

      stop 
      end



