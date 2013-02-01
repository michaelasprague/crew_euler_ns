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

      parameter( nx = 125)
      parameter( ny = 125)
      parameter( nb = 3 )

! solution at time t = 0
      double precision u0  (nx, ny)
      double precision v0  (nx, ny)
      double precision rho0(nx, ny)
      double precision p0  (nx, ny)

! solution at time t 
      double precision u  (nx, ny)
      double precision v  (nx, ny)
      double precision rho(nx, ny)
      double precision p  (nx, ny)

! error

      double precision err_u  (nx, ny)
      double precision err_v  (nx, ny)
      double precision err_p  (nx, ny)
      double precision err_rho  (nx, ny)

! vortex center (true sol)
      double precision vortex_center0(2) !initial
      double precision vortex_center(2)


! solution at time t + tau (acoustic solution)
      double precision ua  (nx, ny)
      double precision va  (nx, ny)
      double precision rhoa(nx, ny)
      double precision pa (nx, ny)

! solution at time t + tau + dtau(acoustic solution)
      double precision ua_new  (nx, ny)
      double precision va_new  (nx, ny)
      double precision rhoa_new(nx, ny)
      double precision pa_new  (nx, ny)

! rhs advective terms
      double precision ru   (nx, ny)
      double precision rv   (nx, ny)
      double precision ruy  (nx, ny)
      double precision rvx  (nx, ny)
      double precision rrhox(nx, ny)
      double precision rrhoy(nx, ny)

! grid stuff
      double precision x  (nx)
      double precision y  (ny)
      double precision xc (nx)
      double precision yc (ny)

! index arrays for simplifying periodic bc implementation
      integer ix ((-nb+1):(nx+nb))
      integer iy ((-nb+1):(ny+nb))

      integer iacoustic(3)
      double precision  dtau(3)

      character*4 filename
      external filename


c======================================================================
c                    RUN PARAMETERS
c======================================================================

      xmin = -10.d0
      xmax = 10.d0 
      ymin = -10.d0
      ymax = 10.d0

      dx   = (xmax-xmin) / float(nx)
      dy   = (ymax-ymin) / float(ny)

      courant  =  0.50d0

      tmax = 10.d0
  
      courant_max = 1.26d0 / sqrt(3.d0) 

      if (courant .gt. courant_max) then
         write(*,*) 'specified courant number too big, resetting' 
         courant = courant_max
         write(*,*) 'courant = ',courant
      endif 

      dt_vtk_out = 1.0d0

      dt_vtk = 0.d0

      i_vtk_out = 0

! initial vortex position
      vortex_center0(1) = -5.d0
      vortex_center0(2) = -5.d0
      

c======================================================================

c----------------------------------------------------------------------
C make the ix,ixu,iy,iyv index arrays; this will simplify
C implementation of the periodic BCs
c     il = 0
c     ir = nx+1
c     do i = 1, nb
c       ix(il) = nx - (i-1)
c       ix(ir) = i 
c       il = il - 1
c       ir = ir + 1
c     enddo
c     do i = 1, nx
c       ix(i) = i
c     enddo
c     il = 0
c     iu = ny+1
c     do i = 1, nb
c       iy(il) = ny - (i-1)
c       iy(iu) = i 
c       il = il - 1
c       iu = iu + 1
c     enddo
c     do i = 1, ny
c       iy(i) = i
c     enddo

      call generate_ix_iy(ix,iy,nb,nx,ny)


c      do i = (-nb+1), (nx+nb)
c         write(*,*) i, ix(i), iy(i)
c      enddo
c      stop

c----------------------------------------------------------------------
C Grid generation
      do i = 1, nx
        x(i)  = xmin + float(i-1)*dx
        xc(i) = x(i) + 0.5d0*dx   ! cell center
      enddo
      do i = 1, ny
        y(i)  = ymin + float(i-1)*dy
        yc(i) = y(i) + 0.5d0*dy   ! cell center
      enddo
c----------------------------------------------------------------------
c initialization

      umean = 1.
      vmean = 1.
      rhomean = 1.

      call init_euler(u0,v0,rho0,x,xc,y,yc,ix,iy,nx,ny,nb,
     &                vortex_center0,umean,vmean,rhomean)

      call update_pressure(p0,rho0,0.,nx,ny)

      call max_speeds_euler(u0, v0, p0, rho0,
     &                      umax, sound_max, nx, ny)

      dt    =  courant * dx / umax

      itmax =  1 + tmax / dt 

      dt = tmax / float(itmax)

      write(*,*) 'dt = ',dt
      write(*,*) 'itmax = ',itmax
 
      write(*,*) 'advection speed', umax
      write(*,*) 'sound speed ', sound_max

      dt_acoustic = dx / (2.d0 * sound_max)
      write(*,*) ' dt / dt_acoustic', dt/dt_acoustic

c     number of acoustic substeps ; must be even 
      ns   = int(dt / dt_acoustic)

      ns_min = 2 

      if (ns .lt. ns_min) then
        ns = ns_min
        write(*,*) 'dt / dt_acoustic too small, setting ns = 2',ns
      else
        write(*,*) 'acoustics substeps; ns = ',ns
      endif

      call copy(u0,   u,  nx,ny)
      call copy(v0,   v,  nx,ny)
      call copy(rho0, rho,nx,ny)
      call copy(p0,   p,  nx,ny)


c      call outputvtk(u0,v0,rho0,p0,xmin,ymin,
c     &     dx,dy,ix,iy,nx,ny,nb,'initial')


c----------------------------------------------------------------------
c number of acoustic substeps for each RK substep
      iacoustic(1) = 1
      iacoustic(2) = ns/2
      iacoustic(3) = ns

      dtau(1) = dt       / 3.d0
      dtau(2) = 0.5d0*dt / float(iacoustic(2))
      dtau(3) = dt       / float(iacoustic(3))

c----------------------------------------------------------------------

      t = 0.d0

c     call r_advect5(ru,ruy,rv,rvx,rrhox,rrhoy,
c    &               u0,v0,rho0,p0,dx,dy,
c    &               ix,iy,nx,ny,nb)

      write(*,*) 'entering main time loop'

      do k = 1, itmax


         call update_euler(t, u, v, rho, p,
     &                     ru, ruy, rv, rvx, rrhox, rrhoy,
     &                     dtau, iacoustic, dx, dy, umax,
     &                     ix, iy, nb, nx, ny, k)

 
        t = t + dt

        write(70,*) t, rho((nx+1)/2,(ny+1)/2)

        write(71,*) t, rho((nx+1)/2,(ny+1)/2)

        if (dt_vtk .gt. dt_vtk_out) then
          dt_vtk = 0.d0
          i_vtk_out = i_vtk_out + 1

c-----true solution at time t-----------------------------------
          vortex_center(1) = vortex_center0(1) + t
          vortex_center(2) = vortex_center0(2) + t
          call init_euler(u0,v0,rho0,x,xc,y,yc,ix,iy,nx,ny,nb,
     &                    vortex_center,umean,vmean,rhomean)
c-----end true solution -----------------------------------------
          call error(u,v,rho,p,u0,v0,rho0,p0,rms_u,rms_rho,
     &               rmsu_rel,
     &               err_max,err_u,err_v,err_rho,err_p,nx,ny)

          option = 1

          call outputvtk(u,v,rho,p,xmin,ymin,dx,dy,
     &       ix,iy,nx,ny,nb,'output'//filename(i_vtk_out),option)

          call outputvtk_err(err_u,err_v,err_rho,p,xmin,ymin,dx,dy,
     &       ix,iy,nx,ny,nb,'output_err'//filename(i_vtk_out))

c----------error-------------------------------
c      call error(u,v,rho,u0,v0,rho0,rms_u,rms_v,rms_rho,
c     &           err_u,err_v,err_rho,nx,ny)
c
c      call outputvtk(err_u,err_v,err_u,p0,xmin,ymin,
c     &     dx,dy,ix,iy,nx,ny,nb,'err_output'//filename(i_vtk_out))
c----------end error----------------------------

        else
          dt_vtk = dt_vtk + dt
        endif

      enddo

      write(*,*) 'finished with main time loop'
      write(*,*) 'tmax = ', t

      call error(u,v,rho,p,u0,v0,rho0,p0,rms_u,rms_rho,
     &           rmsu_rel,err_max,
     &           err_u,err_v,err_rho,err_p,nx,ny)

      write(*,*) 'error: ',rms_u,rms_v,rms_rho

      write(*,*) 'log10(nx), log10(error):'

      write(*,*)  log10(float(nx)), log10(rms_u)

C output -- put everything at cell vertices; this is much averaging
      option =  1
      call outputvtk(u,v,rho,p,xmin,ymin,dx,dy,
     &     ix,iy,nx,ny,nb,'final',option)

      stop
      end

