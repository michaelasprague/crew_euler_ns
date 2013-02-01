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

c code writing started on 05 Nov 2010
c
      program euler2d

      implicit double precision (a-h,o-z)

c ----------------------------------------------------------------------
c  Euler-domain parameters 
c ----------------------------------------------------------------------
c cont0
      parameter(nx_e_quarter = 50) !quarter size of Euler (must be even) 
      parameter(ny_e_quarter = 50) !quarter size of Euler (must be even) 
 
      parameter( nx_e = nx_e_quarter*4+1)  
      parameter( ny_e = ny_e_quarter*4+1)
      parameter( nb_e = 3 )  ! number for 'overlap' for periodic bc's
                             ! rk3 5-6 requires nb_e = 3

c ----------------------------------------------------------------------
c  Navier-Stokes-domain parameters
c ----------------------------------------------------------------------

c nx_n,ny_n tuned to give domain 0.25x0.25 of euler domain
      parameter (nx_n = (nx_e - 1)/4 + 1)  ! for now, should be 
      parameter (ny_n = (ny_e - 1)/4 + 1)
      parameter (mcor_n = 2)   ! number of pressure corrections

      parameter (n_p = 0) !how far projection extends from NS domain
      parameter (nx_p = 2*n_p + nx_n) !projection domain
      parameter (ny_p = 2*n_p + ny_n)

c ----------------------------------------------------------------------
c  Euler-domain arrays
c ----------------------------------------------------------------------

! solution at time t 
      double precision u_e (nx_e, ny_e)
      double precision v_e (nx_e, ny_e)
      double precision r_e (nx_e, ny_e)
      double precision p_e (nx_e, ny_e)

! solution at time t 
      double precision u0_e (nx_e, ny_e)
      double precision v0_e (nx_e, ny_e)
      double precision r0_e (nx_e, ny_e)
      double precision p0_e (nx_e, ny_e)

! error

      double precision err_u_e  (nx_e, ny_e)
      double precision err_v_e  (nx_e, ny_e)
      double precision err_rho_e  (nx_e, ny_e)
      double precision err_p_e  (nx_e, ny_e)

! vortex center (true sol)
      double precision vortex_center0(2) !initial
      double precision vortex_center(2)



! flux terms
      double precision rux (nx_e, ny_e)
      double precision rvy (nx_e, ny_e)
      double precision ruy (nx_e, ny_e)
      double precision rvx (nx_e, ny_e)
      double precision rrx (nx_e, ny_e)
      double precision rry (nx_e, ny_e)

! flux terms
      double precision rux0 (nx_e, ny_e)
      double precision rvy0 (nx_e, ny_e)
      double precision ruy0 (nx_e, ny_e)
      double precision rvx0 (nx_e, ny_e)
      double precision rrx0 (nx_e, ny_e)
      double precision rry0 (nx_e, ny_e)


! grid stuff
      double precision x_e  (nx_e)  ! boundary grid locations
      double precision y_e  (ny_e)
      double precision xc_e (nx_e)  ! center grid locations
      double precision yc_e (ny_e)

! index arrays for simplifying periodic bc implementation
      integer ix ((-nb_e+1):(nx_e+nb_e))
      integer iy ((-nb_e+1):(ny_e+nb_e))

      integer iacoustic(3)
      double precision  dtau(3)

c ----------------------------------------------------------------------
c  Navier-Stokes-domain arrays
c ----------------------------------------------------------------------
      double precision u_n  (nx_n, ny_n)
      double precision v_n  (nx_n, ny_n)
      double precision p_n  (nx_n, ny_n)
      double precision div_n  (nx_n, ny_n)
      double precision r_n  (nx_n, ny_n) ! density in NS domain (all ones)
      double precision energy2d_n  (nx_n, ny_n) !KE density  
      double precision energy2d0_n  (nx_n, ny_n) !KE density  
      double precision err_energy2d_n  (nx_n, ny_n) !KE density error 

c eulers "true" solution on NS domain only (for error tracking)
      double precision u0_n  (nx_n, ny_n)
      double precision v0_n  (nx_n, ny_n)
      double precision p0_n  (nx_n, ny_n)
      double precision r0_n  (nx_n, ny_n)

      double precision ones  (nx_n, ny_n)

c gradient of u and v
      double precision du_dx(nx_n-1,ny_n)
      double precision dv_dx(nx_n-1,ny_n)

      double precision du_dy(nx_n,ny_n-1)
      double precision dv_dy(nx_n,ny_n-1)
     
c error arays
      double precision err_u_n  (nx_n, ny_n)
      double precision err_v_n  (nx_n, ny_n)
      double precision err_p_n  (nx_n, ny_n)

c intermediate arrays used for projection
c      double precision u_proj (nx_n+1, ny_n)
c      double precision v_proj (nx_n, ny_n+1)
      double precision u_proj (nx_p+1, ny_p)
      double precision v_proj (nx_p, ny_p+1)
c      double precision u0_e_proj (nx_n+1, ny_n)
c      double precision v0_e_proj (nx_n, ny_n+1)
      double precision u0_e_proj (nx_p+1, ny_p)
      double precision v0_e_proj (nx_p, ny_p+1)


! velocity boundary conditions (here ns means norht-south)
      double precision u_bc_ew_old(ny_n,2),  u_bc_ns_old(nx_n,2)
      double precision v_bc_ew_old(ny_n,2),  v_bc_ns_old(nx_n,2)

      double precision u_bc_ew(ny_n,2),  u_bc_ns(nx_n,2)
      double precision v_bc_ew(ny_n,2),  v_bc_ns(nx_n,2)

! rho_e, u_e, v_e on NS boundary (calculated in get_navier_bcs)
      double precision re_bc_ew(ny_n,2), re_bc_ns(nx_n,2)
      double precision ue_bc_ew(ny_n,2), ue_bc_ns(nx_n,2)
      double precision ve_bc_ew(ny_n,2), ve_bc_ns(nx_n,2)

! convective bcs
      double precision u_bc_e_conv(ny_n),  u_bc_n_conv(nx_n)
      double precision v_bc_e_conv(ny_n),  v_bc_n_conv(nx_n)

      
      double precision dp_bc_ew(ny_n,2), dp_bc_ns(nx_n,2)
      double precision dp_bc_ew_old(ny_n,2), dp_bc_ns_old(nx_n,2)

      double precision An(nx_n,ny_n) ! A^n times ones 
      double precision A_ew(nx_n-1, ny_n  ) ! interp (a_tilde)^-1 to e-w
      double precision A_ns(nx_n,   ny_n-1) ! interp (a_tilde)^-1 to n-s

! grid stuff
      double precision x_n  (nx_n+1)  ! boundary grid locations
      double precision y_n  (ny_n+1)
      double precision xc_n (nx_n)  ! center grid locations
      double precision yc_n (ny_n)

      double precision tol_cg   (mcor_n)

      double precision b_old( 2*nx_n + 2*ny_n - 4, 2)
      double precision b    ( 2*nx_n + 2*ny_n - 4, 2)

c vtk 
      integer ns_vtk_output ! calc mean (2)or subtract given mean(1)


c ----------------------------------------------------------------------
c  shared
c ----------------------------------------------------------------------
      character*4 filename
      external filename

c     double precision lag_mult(nx_n, ny_n)
      double precision lag_mult(nx_p, ny_p)

c     gradient of lagrange multiplier
      double precision lag_mult_dx(nx_p+1, ny_p)
      double precision lag_mult_dy(nx_p, ny_p+1)

      logical zero_dpdn
      logical forward
      logical two_way
      logical reverse_projection
      logical cond  ! preconditioned CG: on/off
      logical project
      logical neum ! neuman bcs
      logical relative_tol  ! relative(to rhs max norm) tol for press cg
      logical conv    ! convective bcs 
      logical mom     ! momentum field  in proj scheme
      logical matvtk  ! do matlab output
      logical tang    ! multiply tangential bc's by rho_e
      logical ke_flux ! part bound coupling option 
      logical sq      ! part bound coupling option 
      logical mass    ! part bound coupling option 
      logical ke_add  ! run add_ke routine 
 

      open(unit=50,file='params.dat',status='unknown')
      open(unit=51,file='divergence.dat',status='unknown')
      open(unit=52,file='rmsu_n.dat',status='unknown')
      open(unit=53,file='rmsu_e.dat',status='unknown')
      open(unit=54,file='energy.dat',status='unknown')
      open(unit=55,file='b_flux.dat',status='unknown')


c ----------------------------------------------------------------------
c  Constants
c ----------------------------------------------------------------------
       pi = acos(-1.d0)
       zero = 0.
c ----------------------------------------------------------------------
c  User-defined properties for both Euler and NS domains
c ----------------------------------------------------------------------
C c1   CONTROL TOWER 
C c1   control tower
C c1   Control Tower
c    n_p = 0, above for one-way? 

      tmax = 10.d0

      project = .false.

      two_way = .false.

      reverse_projection = .false.

      zero_dpdn = .true.

      cond = .false.  ! preconditioned CG: on/off

      neum = .true. ! Neumann bcs for NS domain (North and East)
c                     make sure n_p =0 above for Neumann

      conv = .false. ! convective bcs for NS domain ( N and E  ) TRASH!

      mom = .true. ! if using mom field instead of velocity field
                   ! its  actually a mass flux not momentum

      tang = .false. ! true for multiplying tangential bc's by rho_e
 
      ke_flux = .false. ! option for neum

      sq  = .true. ! option for neum (currently used)

      mass  = .false. ! option for neum

      ke_add  = .false. ! option for  proj

      matvtk = .false.  ! output matlab format data

c ----------------------------------------------------------------------
c  Euler-domain user-defined quantities
c ----------------------------------------------------------------------
      xmin_e = -10.d0
      xmax_e = 10.d0 
      ymin_e = -10.d0
      ymax_e = 10.d0

      courant_e =  0.5d0

! minimum number of acoustic substeps; must be even
      ns_min = 2 

! initial vortex position
      vortex_center0(1) = -5.d0
      vortex_center0(2) = -5.d0

c ----------------------------------------------------------------------
c  Euler-domain derived quantities
c ----------------------------------------------------------------------
      dx_e = (xmax_e-xmin_e) / float(nx_e)
      dy_e = (ymax_e-ymin_e) / float(ny_e)
  
      courant_max = 1.26d0 / sqrt(3.d0) 

      if (courant_e .gt. courant_max) then
        write(*,*) 'specified courant number too big, resetting' 
        courant_e = courant_max
        write(*,*) 'courant_e = ',courant_e
      endif 

      write(*,*) '---------------------------------'
      write(*,*) 'Euler Domain:'
      write(*,*) 'xmin_e,ymin_e: ',xmin_e,ymin_e
      write(*,*) 'xmax_e,ymax_e: ',xmax_e,ymax_e
      write(*,*) 'dx_e,dy_e: ',dx_e,dy_e

c ----------------------------------------------------------------------
c  Navier-Stokes-domain user-defined quantities
c ----------------------------------------------------------------------
C c2
       
c      viscosity_n = 0.025d0  ! Re = 40
c      viscosity_n = 4.5d-2  ! for grid 134 
c      viscosity_n = 0.09d0  !for grid 68 
      viscosity_n = 0.01d0  ! 
c      viscosity_n = 0.2d0  ! 
c      viscosity_n = 0.4d0  ! 

      density_n   = 1.d0

      relative_tol = .true. !for pressure solve
      tol_cg(1) = 1.d-8
      tol_cg(2) = 1.d-8

      tol_bicg = 1.d-8
 
      tol_bicg_p = 1.d-8 ! tol for unsymetric press (not used)

      tol_cg_proj = 1.d-8 ! tol foe cg_solve in projection scheme
                           ! does not use relative tol


c ----------------------------------------------------------------------
c  Navier-Stokes-domain derived quantities
c ----------------------------------------------------------------------

c the following assumes matching grids; need to change

c      nx_start = (nx_e - 1) * 3/ (4 * 2)
c      ny_start = (ny_e - 1) * 3/ (4 * 2)

      nx_start = (nx_e-1) / 2 - (nx_n-1)/2 + 1
      ny_start = (ny_e-1) / 2 - (ny_n-1)/2 + 1

      nx_p_start = nx_start - n_p
      ny_p_start = ny_start - n_p

      write (*,*) 'nx_e',nx_e
      write (*,*) 'nx_n',nx_n
      write (*,*) 'nx_p',nx_p
      write (*,*) 'nx_start',nx_start
      write (*,*) 'ny_start',ny_start
      write (*,*) 'nx_p_start',nx_p_start
      write (*,*) 'ny_p_start',ny_p_start
      write (*,*) '(nx_n-1)/2', (nx_n-1)/2
      write (*,*) 'float (nx_n-1)/2', float((nx_n-1)/2)
c      stop

c     xmin_n = (nx_start - 1) * dx_e
c     ymin_n = (ny_start - 1) * dy_e
c     xmax_n = xmin_n + nx_n*dx_e
c     ymax_n = ymin_n + ny_n*dy_e

      xmin_n = -dx_e * (float((nx_n-1)/2) + 0.5d0)
      ymin_n = -dy_e * (float((ny_n-1)/2) + 0.5d0)
      xmax_n =  dx_e * (float((nx_n-1)/2) + 0.5d0)
      ymax_n =  dy_e * (float((ny_n-1)/2) + 0.5d0)

c assuming same grid spacing for both euler and ns; will change later

      dx_n = (xmax_n - xmin_n) / float(nx_n)
      dy_n = (ymax_n - ymin_n) / float(ny_n)

      write(*,*) '---------------------------------'
      write(*,*) 'Navier-Stokes Domain:'
      write(*,*) 'xmin_n,ymin_n: ',xmin_n,ymin_n
      write(*,*) 'xmax_n,ymax_n: ',xmax_n,ymax_n
      write(*,*) 'dx_n,dy_n: ',dx_n,dy_n
      write(*,*) 'nx_start, ny_start: ',nx_start,ny_start
      write(*,*) '---------------------------------'
    

c ----------------------------------------------------------------------
c                    RUN PARAMETERS
c ----------------------------------------------------------------------
c  c3   vtk control
     
      dt_vtk_out =0.0d0

      dt_vtk = 0.d0
      i_vtk_out = 0

      
      ns_vtk_option = 1 ! use given umean and vmean
c      ns_vtk_option = 2 ! calc mean inside 


c mean velocity for outputvtk_ns (subtracts it)
c and for init_euler

      umean = 1.d0
      vmean = 1.d0
      rhomean = 1.d0
      pmean = 1.d0 
 
c this is how it looks down there in vtk output
c       if (dt_vtk .ge. dt_vtk_out) then
c          dt_vtk = 0.d0
c          i_vtk_out = i_vtk_out + 1
c          call outputvtk(....i_vtk_out....
c        else
c          dt_vtk = dt_vtk + dt
c        endif

c-----------------------------------------------------------------------


c======================================================================

! create ix & iy arrays, which are needed for periodic boundaries
      call generate_ix_iy(ix, iy, nb_e, nx_e, ny_e)

c----------------------------------------------------------------------
C Grid generation: Euler
      do i = 1, nx_e
        x_e(i)  = xmin_e + float(i-1)*dx_e
        xc_e(i) = x_e(i) + 0.5d0*dx_e   ! cell center
      enddo
      do i = 1, ny_e
        y_e(i)  = ymin_e + float(i-1)*dy_e
        yc_e(i) = y_e(i) + 0.5d0*dy_e   ! cell center
      enddo
c----------------------------------------------------------------------
C Grid generation: NS
      do i = 1, nx_n
        x_n(i)  = xmin_n + float(i-1)*dx_n
        xc_n(i) = x_n(i) + 0.5*dx_n   ! cell center
      enddo
      x_n(nx_n+1)  = xmax_n
      do i = 1, ny_n
        y_n(i)  = ymin_n + float(i-1)*dy_n
        yc_n(i) = y_n(i) + 0.5*dy_n   ! cell center
      enddo
      y_n(ny_n+1)  = ymax_n

c------------debug--------------
c      do i = 1,nx_e
c         write(*,*) 'x_e(',i,')',x_e(i)
c      enddo
c
c      do i = 1,nx_e
c         write(*,*) 'xc_e(',i,')',xc_e(i)
c      enddo
c-----------end debug--------------

       write(*,*) 'NS domain limits'
       write(*,*) 'x_n(1), x_n(nx_n+1)',x_n(1), x_n(nx_n+1)
       write(*,*) 'y_n(1), y_n(ny_n+1)',y_n(1), y_n(ny_n+1)
c       stop

c----------------------------------------------------------------------
c initialize arrays
c----------------------------------------------------------------------

      call init_euler(u_e, v_e, r_e, x_e, xc_e, y_e, yc_e,
     &                      ix, iy, nx_e, ny_e, nb_e, vortex_center0,
     &                       umean, vmean, rhomean)

! putting the same thing in to u0 arrays; used later for error checking
c      call init_euler(u0_e, v0_e, r0_e, x_e, xc_e, y_e, yc_e,
c     &                      ix, iy, nx_e, ny_e, nb_e, vortex_center0,
c     &                       umean, vmean, rhomean)

      
      call update_pressure(p_e, r_e, zero, nx_e, ny_e)
c      call update_pressure(p0_e, r0_e, 0., nx_e, ny_e)


      call max_speeds_euler(u_e, v_e, p_e, r_e,
     &                      u_max, sound_max, nx_e, ny_e)

c calculate dt based on grid, max velocity, and user-spec courant
      dt    =  courant_e * dx_e / u_max


c calculate dt based on stability requirement
      dt = min(2.*viscosity_n / u_max**2 , 0.5*dx_e**2 / viscosity_n)
      dt = courant_e * dt

      itmax =  1 + int(tmax / dt)

      dt = tmax / float(itmax)

      dt_acoustic = dx_e / (2.d0 * sound_max)

c     number of acoustic substeps ; must be even 
      ns   = int(dt / dt_acoustic)
 
      if (mod(ns,2) .ne. 0) ns = ns + 1      

      if (ns .lt. ns_min) then
        ns = ns_min
        write(*,*) 'dt / dt_acoustic too small, setting ns = ',ns
      endif

c      write(*,*) 'dt=', dt, 'itmax=',itmax
c      stop

c outputvtk subtracting mean:
c 1 for mean = 1; 
c 2 for mean = 0;
c      call outputvtk(u_e,v_e,r_e,p_e,xmin_e, ymin_e,dx_e,dy_e,
c     &       ix,iy,nx_e,ny_e,nb_e,'output'//filename(0),1)

c----------------------------------------------------------------------
c number of acoustic substeps for each RK substep
      iacoustic(1) = 1
      iacoustic(2) = ns/2
      iacoustic(3) = ns

      dtau(1) = dt       / 3.d0
      dtau(2) = 0.5d0*dt / float(iacoustic(2))
      dtau(3) = dt       / float(iacoustic(3))

c----------------------------------------------------------------------
c write out interesting run parameters
      write(*,*) 'dt = ',dt
      write(*,*) 'itmax = ',itmax
      write(*,*) 'acoustics substeps; ns = ',ns
      write(*,*) ' dt / dt_acoustic', dt/dt_acoustic
 
      write(*,*) 'nx_e, ny_e', nx_e, ny_e
      write(*,*) 'nx_n, ny_n', nx_n, ny_n

c----------------------------------------------------------------------


c----------write run parameters into file params.dat--------------------


       write(50,*) 'nx_e,ny_e',nx_e,ny_e
       write(50,*) 'nx_n,ny_n',nx_n,ny_n
       write(50,*) 'xmin_n,ymin_n', xmin_n, ymin_n
       write(50,*) 'xmax_n,ymax_n', xmax_n, ymax_n
       write(50,*) 'dt',dt
       write(50,*) 'dx_n,dy_n',dx_n,dy_n
       write(50,*) 'dx_e,dy_e',dx_e,dy_e
       write(50,*) 'n_p',n_p
       write(50,*) 'dt_vtk_out',dt_vtk_out
       write(50,*) 'viscosity_n',viscosity_n
       write(50,*) 'u_max_euler',u_max
       write(50,*) 'Peclet',(u_max*dx_e)/viscosity_n
       write(50,*) 'Peclet * CFL',(u_max*dx_e)/viscosity_n
       write(50,*) 'project',project
       write(50,*) 'two_way',two_way
       write(50,*) 'zero_dpdn',zero_dpdn
       write(50,*) 'neum',neum
       write(50,*) 'reverse_projection',reverse_projection
       write(50,*) 'mom',mom
       write(50,*) 'tang',tang
       write(50,*) 'sq',sq
       write(50,*) 'mass',mass
       write(50,*) 'ke_flux',ke_flux
       write(50,*) 'ke_add',ke_add
       call flush(50)


c--------------------end write parameters----------------------


      t = 0.d0

c     itmax = 2

      do k = 1, itmax
c5
        write(*,*) '**************************************************'
        write(*,*) '        t = ', t,'i_t=',k,'i_vtk_out',i_vtk_out
        write(*,*) '**************************************************'
        call update_euler(t, u_e, v_e, r_e, p_e,
     &                    rux, ruy, rvy, rvx, rrx, rry,
     &                    dtau, iacoustic, dx_e, dy_e, u_max,
     &                    ix, iy, nb_e, nx_e, ny_e, k)

c-------------euler--"true"--solution------------------------------
c        call update_euler(t, u0_e, v0_e, r0_e, p0_e,
c     &                    rux0, ruy0, rvy0, rvx0, rrx0, rry0,
c     &                    dtau, iacoustic, dx_e, dy_e, u_max,
c     &                    ix, iy, nb_e, nx_e, ny_e, k)
c-------------end--euler--"true"--solution------------------------------

c        call grab_bound_vels(u_e, v_e, u_proj, v_proj, 
c     &                       nx_e, ny_e, nx_n, ny_n,
c     &                       nx_start, ny_start)

        call grab_bound_vels(u_e, v_e, r_e, u_proj, v_proj,
     &                       nx_e, ny_e, nx_p, ny_p,
     &                       nx_p_start, ny_p_start, mom)

c----debug------------
c        do i =1,nx_n
c        write(*,*) 'before', u_proj(1,i)-1.d0 + v_proj(i,1)-1.d0
cc         write(*,*) 'befor r', r_e(nx_p_start+i,ny_p_start) -
cc     &                          r_e(nx_p_start,i+ny_p_start)
c        enddo
c         write(*,*) 'vs,uw,rs,ew', v_e(nx_start+2,ny_start),
c     &                             u_e(nx_start,ny_start+2),
c     &                             r_e(nx_start+2,ny_start),
c     &                             r_e(nx_start,ny_start+2),
c     &                             v_proj(3,1),u_proj(1,3)
         
c---end-debug----------

       if (project) then
        forward = .true. 
c        call projection(u_proj, v_proj, lag_mult, dx_n, dy_n, 
c     &                  nx_n + 1, ny_n + 1, forward)


        call projection(u_proj, v_proj, lag_mult, dx_n, dy_n,
     &                  nx_p, ny_p, forward, mom, density_n,
     &                  lag_mult_dx, lag_mult_dy, tol_cg_proj)

c--------start - debug--------------------
c--------end   - debug--------------------


c        call get_max_div(div_max,bd_flux, u_proj, v_proj, dx_n, dy_n, 
c     &                   nx_n + 1, ny_n + 1)

        call get_max_div(div_max,bd_flux, u_proj, v_proj, dx_n, dy_n,
     &                   nx_p + 1, ny_p + 1)

        write(51,*) t, div_max, bd_flux

        write(*,*) 'div after proj: ',div_max
c           write(*,*) 'k = ', k, 't = ', t
        if (div_max .gt. 1.d-2) then
           write(*,*) 'k = ', k, 't = ', t
           stop 'div too great'
        endif
        
       endif !projection

c----------------------------------------------------------------------
        if (k .eq. 1) then

          write(*,*) 'initializing NS domain'

          call init_navier_stokes(u_n, v_n, p_n, u_proj, v_proj, 
     &                            dx_n, dy_n, nx_n, ny_n, n_p)

c          call outputvtk_ns(u_n, v_n, p_n, p_n, xmin_n, ymin_n,
c     &                      dx_n, dy_n, nx_n, ny_n,
c     &                      'output_ns'//filename(0),
c     &                       umean,vmean,ns_vtk_option)


c if (tang) then tangential bc's are multiplied by rho
c normal velocities are multiplied by rho in
c grab_bound_vels routine if (mom) 


          call get_navier_bcs(u_e, v_e, r_e, u_proj, v_proj,
     &               u_bc_ew_old, u_bc_ns_old, v_bc_ew_old, v_bc_ns_old,
     &               ue_bc_ew, ue_bc_ns, ve_bc_ew, ve_bc_ns,
     &               re_bc_ns, re_bc_ew,
     &               dp_bc_ew_old, dp_bc_ns_old,
     &               nx_start, ny_start,
     &               nx_e, ny_e, nx_n, ny_n, n_p, tang, density_n,
     &               neum, ke_flux,sq, mass,project,
     &               lag_mult, dx_n, dy_n)



c----------add KE---------------------------------------
          if (ke_add) then
            call add_ke3( u_bc_ew_old, u_bc_ns_old,
     &                    v_bc_ew_old, v_bc_ns_old,
     &                   ue_bc_ew, ue_bc_ns, ve_bc_ew, ve_bc_ns,
     &                   re_bc_ns, re_bc_ew,
     &                   nx_n,ny_n,dx_n,dy_n,density_n, alpha)
          endif
c---------end-add KE-----------------------------------       


          call  outputvtk_bound(u_bc_ew_old, u_bc_ns_old,
     &                 v_bc_ew_old, v_bc_ns_old,
     &                 xmin_n, ymin_n, xmax_n, ymax_n,
     &                 umean, vmean,
     &                 dx_n,dy_n,nx_n,ny_n, 
     &                 'output_bound'//filename(0))

          if (conv) then
            call get_conv_bcs(u_bc_ew_old, u_bc_ns_old,
     &                        v_bc_ew_old, v_bc_ns_old,
     &                        u_bc_e_conv, u_bc_n_conv,
     &                        v_bc_e_conv, v_bc_n_conv,
     &                        u_n,v_n,
     &                        dt, nx_n, ny_n, dx_n, dy_n, k)
          endif
  
          call calc_b(b_old,
     &          u_bc_ew_old, u_bc_ns_old,
     &          v_bc_ew_old, v_bc_ns_old,
     &          u_bc_e_conv, u_bc_n_conv,
     &          v_bc_e_conv, v_bc_n_conv,
     &          nx_n, ny_n, dx_n, dy_n, viscosity_n, neum, conv)

          write(*,*) '1',dy_n
 
          call get_pbc_navier_stokes(u_n,v_n,
     &            u_bc_ew_old, u_bc_ns_old, v_bc_ew_old, v_bc_ns_old,
     &            u_bc_ew_old, u_bc_ns_old, v_bc_ew_old, v_bc_ns_old,
     &            dp_bc_ew_old, dp_bc_ns_old, 
     &            dt, nx_n, ny_n, 
     &            dx_n, dy_n, density_n, viscosity_n,zero_dpdn)

        endif
c----------------------------------------------------------------------



c if (tang) then tangential bc's are multiplied by rho
c normal velocities are multiplied by rho in
c grab_bound_vels routine if (mom) 
        call get_navier_bcs(u_e, v_e, r_e, u_proj, v_proj,
     &       u_bc_ew, u_bc_ns, v_bc_ew, v_bc_ns,
     &       ue_bc_ew, ue_bc_ns, ve_bc_ew, ve_bc_ns,
     &       re_bc_ns, re_bc_ew,
     &       dp_bc_ew, dp_bc_ns,
     &       nx_start, ny_start,
     &       nx_e, ny_e, nx_n, ny_n, n_p, tang, density_n,
     &       neum, ke_flux, sq, mass, project,
     &       lag_mult, dx_n, dy_n)

c7---------debug: write bc's to a file------------      
c      if (k .eq. 304) then
c      do i = 1, ny_n
c        write(55,*) i,u_bc_ew(i,1),u_bc_ew(i,2)
c      enddo
c      endif
c---------end- debug: write bc's to a file------------      


c--------start-add KE---------------------------------------
        if (ke_add) then
          call add_ke3( u_bc_ew, u_bc_ns, v_bc_ew, v_bc_ns,
     &                 ue_bc_ew, ue_bc_ns, ve_bc_ew, ve_bc_ns,
     &                 re_bc_ns, re_bc_ew,
     &                 nx_n,ny_n,dx_n,dy_n,density_n,alpha)
        endif
c---------end-add KE-----------------------------------       


c-------start - debug--------------------
c-------end-debug------------------
        if (conv) then
          call get_conv_bcs(u_bc_ew, u_bc_ns,
     &                      v_bc_ew, v_bc_ns,
     &                      u_bc_e_conv, u_bc_n_conv,
     &                      v_bc_e_conv, v_bc_n_conv,
     &                      u_n,v_n,
     &                      dt, nx_n, ny_n, dx_n, dy_n, k)

        endif

        write(*,*) 'START UPDATE NS'


        call update_navier_stokes(t, u_n, v_n, p_n, div_n,
     &          u_bc_ew_old,  u_bc_ns_old,  v_bc_ew_old,  v_bc_ns_old,
     &          u_bc_ew,      u_bc_ns,      v_bc_ew,      v_bc_ns,
     &          u_bc_e_conv, u_bc_n_conv, v_bc_e_conv, v_bc_n_conv,
     &          dp_bc_ew,     dp_bc_ns, 
     &          b_old, b,
     &          An, A_ew, A_ns,
     &          tol_bicg, tol_cg,
     &          dx_n, dy_n, dt, viscosity_n, density_n, mcor_n,
     &          nx_n, ny_n, k, zero_dpdn, cond, neum, relative_tol,conv)


c8------debug: calc boundary flux------------------------
c        b_flux = 0.d0
c        do i = 1,ny_n
c           b_flux = b_flux - u_bc_ew(i,1)*dy_n + u_bc_ew(i,2)*dy_n
c        enddo
c        do i = 1,nx_n
c           b_flux = b_flux - v_bc_ns(i,1)*dx_n + v_bc_ns(i,2)*dx_n
c        enddo
c        if (abs(b_flux) .gt. 10.*tol_cg_proj) then
c           write(*,*) 'b_flux exceeds limit:', b_flux
c           stop ' b_flux exceeds limit'
c        endif
        !write(55,*) t,abs(b_flux)
        !call flush(55)

c------end debug: calc boundary flux------------------------

c------debug: calc divergence after NS step---------------


c------debug: calc divergence after NS step---------------

        write(*,*) 'END UPDATE NS'

        if (two_way) then
          write(*,*) 'START TWO-WAY'

C should the following be rhie-chow velocities?
          do j = 1, ny_n
            do i = 2, nx_n
c              u_proj(i,j) = 0.5d0 * (u_n(i-1,j) + u_n(i,j))
              u_proj(i+n_p,j+n_p) = 0.5d0 * (u_n(i-1,j) + u_n(i,j))
            enddo
          enddo
          do i = 1, nx_n
            do j = 2, ny_n
c              v_proj(i,j) = 0.5d0 * (v_n(i,j-1) + v_n(i,j))
              v_proj(i+n_p,j+n_p) = 0.5d0 * (v_n(i,j-1) + v_n(i,j))
            enddo
          enddo

          if (reverse_projection) then
            forward = .false. 
            call projection(u_proj, v_proj, lag_mult, dx_n, dy_n, 
     &                      nx_p, ny_p, forward, mom, density_n,
     &                      lag_mult_dx, lag_mult_dy,tol_cg_proj)
  
          endif

          call put_bound_vels(u_e, v_e,r_e, u_proj, v_proj, 
     &                        nx_e, ny_e, nx_p, ny_p,
     &                        nx_p_start, ny_p_start, mom) 


          write(*,*) 'END TWO-WAY'

        endif

 
        t = t + dt

        write(70,*) t, r_e((nx_e+1)/2,(ny_e+1)/2)
c------------------------------------------------------------------------
c                 START VTK
c------------------------------------------------------------------------
c c4
        if (dt_vtk .ge. dt_vtk_out) then
          dt_vtk = 0.d0
          i_vtk_out = i_vtk_out + 1
        
          write(*,*) 'START VTK'

c EULER
c outputvtk subtracting mean:
c 1 for mean = 1; 
c 2 for mean = 0;

c         call outputvtk(u_e,v_e,r_e,p_e,xmin_e, ymin_e,dx_e,dy_e,
c     &                  ix,iy,nx_e,ny_e,nb_e,
c     &                 'output_e'//filename(i_vtk_out),1)

c         call outputvtk(u_e,v_e,r_e,p_e,xmin_e, ymin_e,dx_e,dy_e,
c     &                  ix,iy,nx_e,ny_e,nb_e,
c     &                 'output_e_nomean'//filename(i_vtk_out),2)
c
c         call outputvtk(u0_e,v0_e,r0_e,p0_e,xmin_e, ymin_e,dx_e,dy_e,
c     &                  ix,iy,nx_e,ny_e,nb_e,
c     &                 'output_e_true'//filename(i_vtk_out),.false.)

c         call error(u_e,v_e,r_e,p_e,u0_e,v0_e,r0_e,p0_e,rmsu_e,rmsr_e,
c     &           rmsu_rel_e,err_max_e,err_u_e,err_v_e,err_rho_e,err_p_e,
c     &            nx_e,ny_e)
c
c         call outputvtk_err(err_u_e,err_v_e,err_rho_e,err_p_e,
c     &                      xmin_e,ymin_e,
c     &                      dx_e,dy_e,ix,iy,nx_e,ny_e,nb_e,
c     &                      'output_e_err'//filename(i_vtk_out))

c        write(53,*) t,rmsu_e,rmsu_rel_e,err_max_e
       

c NAVIER STOKES 


c----------true solution on NS domain only--------------------------
c        call grab_bound_vels(u0_e, v0_e, r0_e, u0_e_proj, v0_e_proj,
c     &                       nx_e, ny_e, nx_p, ny_p,
c     &                       nx_p_start, ny_p_start,.false.)
c        
c        call init_navier_stokes(u0_n, v0_n, p0_n,
c     &                          u0_e_proj, v0_e_proj,
c     &                          dx_n, dy_n, nx_n, ny_n, n_p) !  


          call copy_e_to_n(u_e, v_e, r_e,
     &                    u0_n, v0_n, r0_n,
     &                    nx_e, ny_e, nx_n, ny_n,
     &                    nx_start, ny_start)


c----------end true solution on NS domail only---------------------

c---------------kinetic energy-----------
c energy routines 

          call initialize(density_n,r_n,nx_n,ny_n)

          call calc_energy_n(energy_n, energy2d_n, u_n, v_n, r_n,
     &                     dx_n, dy_n, nx_n, ny_n)
          call calc_energy_n(energy_n0, energy2d0_n, u0_n, v0_n, r0_n,
     &                     dx_n, dy_n, nx_n, ny_n)
c         call calc_energy_e(energy_e,u_e,v_e,r_e,dx_e,dy_e,nx_e,ny_e)

          call error_ns_ke(energy2d_n, energy2d0_n, err_energy2d_n,
     &                 rmse,rmse_rel,
     &                 nx_n, ny_n)


          write(*,*)'energy', i_vtk_out, t, energy_n, energy_n0, alpha


          write(54,*) i_vtk_out, t, energy_n, energy_n0, alpha
          call flush(54)

c--------------end kinetic energy-------------

c NS
c outputvtk_ns: option =1 - subtract given mean
c               option =2 - calc mean inside


c       call outputvtk_ns(u_n, v_n, p_n, energy2d_n,
c    &                      xmin_n, ymin_n, dx_n, dy_n, nx_n, ny_n,
c    &                      'output_ns'//filename(i_vtk_out),
c    &                      umean,vmean,pmean,ns_vtk_option)

          call outputvtk_ns(u_n, v_n, p_n, div_n,
     &                      xmin_n, ymin_n, dx_n, dy_n, nx_n, ny_n,
     &                      'output_ns'//filename(i_vtk_out),
     &                      umean,vmean,pmean,ns_vtk_option)


          call  outputvtk_bound(u_bc_ew, u_bc_ns,
     &                 v_bc_ew, v_bc_ns,
     &                 xmin_n, ymin_n, xmax_n, ymax_n,
     &                 umean, vmean,
     &                 dx_n,dy_n,nx_n,ny_n,
     &                 'output_bound'//filename(i_vtk_out))


          call initialize(density_n,r_n,nx_n,ny_n)

c if (mom) then mass flux is used for error calc

          call error_ns(u_n, v_n, r_n, u0_n, v0_n, r0_n,err_u_n,err_v_n,
     &                 rmsu_n, rmsu_rel_n,err_max_n,
     &                nx_n, ny_n, mom)

          write(52,*) t,rmsu_n,rmsu_rel_n,err_max_n
         

c-------------------------------debug-----------------------

c------------end- grad u and grad v------------


c--------------------------------------------------------------------------------
c Lagrange multiplier and its gradient
c
c      call matlabvtk(lag_mult,nx_p,ny_p,
c     &              'lag_mult'//filename(i_vtk_out))
c 
c      call matlabvtk(lag_mult_dx,nx_p+1,ny_p,
c     &              'lag_mult_dx'//filename(i_vtk_out))
c
c      call matlabvtk(lag_mult_dy,nx_p,ny_p+1,
c     &              'lag_mult_dy'//filename(i_vtk_out))
c-------------------------------------------------------------------------------

c--------------------------------------end-debug-----------------------------

c MATLAB output of velocities -- don't trust Paraview (bs'd me multiple times)

          if (matvtk) then

            write(*,*) 'START MATLAB VTK'


c ENERGY

            call matlabvtk(energy2d_n, nx_n, ny_n,
     &              'energy2d_n'//filename(i_vtk_out))

c NS
            call matlabvtk(u_n,nx_n,ny_n,
     &              'u_vel_n'//filename(i_vtk_out))
            call matlabvtk(v_n,nx_n,ny_n,
     &               'v_vel_n'//filename(i_vtk_out))

c EULER

            call matlabvtk(u_e,nx_e,ny_e,
     &               'u_vel_e'//filename(i_vtk_out))
            call matlabvtk(v_e,nx_e,ny_e,
     &               'v_vel_e'//filename(i_vtk_out))

c NS Boundary 
c           call matlabvtk_bound(u_bc_ew, u_bc_ns,
c    &                      v_bc_ew, v_bc_ns,
c    &                      nx_n, ny_n, 'ns_bound'//filename(i_vtk_out))



            write(*,*) 'END MATLAB VTK'

          endif ! matlabvtk

c end MATLAB VTK
c------------------------------------------------------------------------------  
 
          write(*,*) 'END VTK'

        else
          dt_vtk = dt_vtk + dt
        endif
c------------------------------------------------------------------------
c               END VTK
c------------------------------------------------------------------------

      enddo !time loop

      write(*,*) 'tmax = ', t
      write(*,*) 'dt = ', dt
      write(*,*) 'i_t = ',k 
      write(*,*) 'nx_e = ',nx_e,'nx_n',nx_n

c      call error(u_e,v_e,r_e,p_e,u0_e,v0_e,r0_e,p0_e,rmsu_e,rmsr_e,
c     &          rmsu_rel_e,err_max_e,err_u_e,err_v_e,err_rho_e,err_p_e,
c     &          nx_e,ny_e)

c     write(*,*) 'error_e: ','rmsu_e',rmsu_e,'rmsu_rel_e',rmsu_rel_e

c      write(*,*) 'log10(nx), log10(error):'
c      write(*,*)  log10(float(nx_e)), log10(rms_u)

C output -- put everything at cell vertices; this is much averaging

c      call outputvtk(u_e,v_e,r_e,p_e,xmin_e, ymin_e,dx_e,dy_e,
c     &               ix,iy,nx_e,ny_e,nb_e,'final')



      stop
      end

