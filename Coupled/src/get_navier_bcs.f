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

      subroutine get_navier_bcs(u_e, v_e, r_e, u_proj, v_proj,
     &                          u_bc_ew, u_bc_ns, v_bc_ew, v_bc_ns,
     &                          ue_bc_ew, ue_bc_ns, ve_bc_ew, ve_bc_ns,
     &                          re_bc_ns, re_bc_ew,
     &                          dp_bc_ew, dp_bc_ns,
     &                          nx_start, ny_start,
     &                          nx_e,ny_e,nx_n,ny_n,n_p,tang,density,
     &                          neum,ke_flux,sq,mass,project,
     &                          p, dx_n, dy_n)

! this takes the projected divergence-free normal velocities along with the
! original euler velocities and fills the boundary condition arrays for the 
! Navier-STokes solver.
!
! This is a little complicated, because we need to interpolate the euler
! velocities to get appropriate tangential velocities, which, of course, do
! not come into play in the projection scheme

c if (mom) then tangential velocities are multiplied by rho


      implicit double precision (a-h,o-z)

      double precision u_e (nx_e, ny_e)
      double precision v_e (nx_e, ny_e)
      double precision r_e (nx_e, ny_e)

      double precision u_proj (nx_n+2*n_p+1, ny_n+2*n_p)
      double precision v_proj (nx_n+2*n_p, ny_n+2*n_p+1)

      double precision  p(nx_n+2*n_p, ny_n+2*n_p) !lagrange mult

      double precision u_bc_ew(ny_n,2), u_bc_ns(nx_n,2)
      double precision v_bc_ew(ny_n,2), v_bc_ns(nx_n,2)
      double precision dp_bc_ew(ny_n,2), dp_bc_ns(nx_n,2)

      double precision re_bc_ew(ny_n,2), re_bc_ns(nx_n,2)
      double precision ue_bc_ew(ny_n,2), ue_bc_ns(nx_n,2)
      double precision ve_bc_ew(ny_n,2), ve_bc_ns(nx_n,2)

c inside arrays
      double precision au_bc_ew(ny_n,2), au_bc_ns(nx_n,2)
      double precision av_bc_ew(ny_n,2), av_bc_ns(nx_n,2)
  
      double precision p_x(nx_n+2*n_p) !interpol lambda horiz  
      double precision p_y(ny_n+2*n_p) !interpol lambda vert
 
      logical tang 
      logical neum 
      logical ke_flux 
      logical sq 
      logical mass 
      logical project 


      nx_p = nx_n + 2*n_p
      ny_p = ny_n + 2*n_p





! south

c--------------interpolate lambd horiz---------------
      if (n_p .ge. 1) then
        do i = 1,nx_p
           p_x(i) = 0.5d0*( p(i,n_p) + p(i,n_p+1))          
         enddo
      endif
 
c--------------end interp lambda--------------------



      do i = 1, nx_n

c--------------copy u_e, v_e, rho_e to NS boundary----------------------

        re_bc_ns(i,1) = 0.5d0*
     &                   ( r_e(nx_start+i-1,ny_start) +
     &                     r_e(nx_start+i-1,ny_start-1) )

        
        vel1 = u_e(nx_start+i-1, ny_start  )
        vel2 = u_e(nx_start+i-1, ny_start-1)
        vel3 = u_e(nx_start+i,   ny_start  )
        vel4 = u_e(nx_start+i,   ny_start-1)

        u_bc_ns(i,1) = 0.25*(vel1+vel2+vel3+vel4)

        
        v_bc_ns(i,1) = v_e(nx_start+i-1,ny_start)


c---copy u_e,v_e for sending into add_ke routine-------
         ue_bc_ns(i,1) = u_bc_ns(i,1)
         ve_bc_ns(i,1) = v_bc_ns(i,1)


c------------------end-copy---------------------------------------------- 
              
        if (tang .and. project) then
          u_bc_ns(i,1) = u_bc_ns(i,1)*re_bc_ns(i,1)
        endif

 
        if (project) then
          v_bc_ns(i,1) = v_proj(i+n_p,1+n_p)
   
          if (n_p .ge. 1) then 
          u_bc_ns(i,1) = u_bc_ns(i,1) -
     &                   (p_x(n_p+2 +i-1) - p_x(n_p +i-1))/(2.d0*dx_n)
          endif

        endif

        if (tang .and. project) then
          u_bc_ns(i,1) = u_bc_ns(i,1)/density
        endif


        if (ke_flux .and. neum) then

           au_bc_ns(i,1) = (re_bc_ns(i,1)/density)*
     &                    ( u_bc_ns(i,1)**2 + v_bc_ns(i,1)**2 )*
     &                      u_bc_ns(i,1)
 
           av_bc_ns(i,1) = (re_bc_ns(i,1)/density)*
     &                    ( u_bc_ns(i,1)**2 + v_bc_ns(i,1)**2 )*
     &                      v_bc_ns(i,1)

           u_bc_ns(i,1) = au_bc_ns(i,1)*
     &                  ( au_bc_ns(i,1)**2 + av_bc_ns(i,1)**2)**
     &                  (-1.d0/3.d0)
  
           v_bc_ns(i,1) = av_bc_ns(i,1)*
     &                  ( au_bc_ns(i,1)**2 + av_bc_ns(i,1)**2)** 
     &                  (-1.d0/3.d0)
  
         endif ! ke_flux and neum

         if (sq .and. neum) then

          u_bc_ns(i,1) = u_bc_ns(i,1)*sqrt(re_bc_ns(i,1)/density)
          v_bc_ns(i,1) = v_bc_ns(i,1)*sqrt(re_bc_ns(i,1)/density)

         endif ! sq and neum

         if (mass .and. neum) then

          u_bc_ns(i,1) = u_bc_ns(i,1)*(re_bc_ns(i,1)/density)
          v_bc_ns(i,1) = v_bc_ns(i,1)*(re_bc_ns(i,1)/density)

         endif ! mass and neum


      enddo ! i=1,nx_n 
 
 

! north

c--------------interpolate lambd horiz---------------
      if (n_p .ge. 1) then
        do i = 1,nx_p
           p_x(i) = 0.5d0*( p(i,n_p+ny_n) + p(i,n_p+ny_n+1))
        enddo
      endif

c--------------end interp lambda--------------------



      do i=1,nx_n
c--------------copy u_e, v_e, rho_e to NS boundary----------------------

        re_bc_ns(i,2) = 0.5d0*
     &                   ( r_e(nx_start+i-1,ny_start+ny_n) +
     &                     r_e(nx_start+i-1,ny_start+ny_n-1) )
 
        vel1 = u_e(nx_start+i-1, ny_start+ny_n  )
        vel2 = u_e(nx_start+i-1, ny_start+ny_n-1)
        vel3 = u_e(nx_start+i,   ny_start+ny_n  )
        vel4 = u_e(nx_start+i,   ny_start+ny_n-1)

        u_bc_ns(i,2) = 0.25*(vel1+vel2+vel3+vel4)
       
        v_bc_ns(i,2) = v_e(nx_start+i-1,ny_start+ny_n)

c---copy u_e,v_e for sending into add_ke routine-------
         ue_bc_ns(i,2) = u_bc_ns(i,2)
         ve_bc_ns(i,2) = v_bc_ns(i,2)

c---------------------end-copy-------------------------------------------- 
       
        if (tang .and. project) then
          u_bc_ns(i,2) = u_bc_ns(i,2)*re_bc_ns(i,2)
        endif

 
        if (project) then 
          v_bc_ns(i,2) = v_proj(i+n_p,ny_n+n_p+1)
          
          
          if (n_p .ge. 1) then
          u_bc_ns(i,2) = u_bc_ns(i,2) -
     &                   (p_x(n_p+2 +i-1) - p_x(n_p +i-1))/(2.d0*dx_n)
          endif


        endif

        if (tang .and. project) then
          u_bc_ns(i,2) = u_bc_ns(i,2)/density
        endif



      enddo ! i=1,nx_n 


! west

c--------------interpolate lambd vertic---------------

      if (n_p .ge. 1) then
        do j = 1,ny_p
           p_y(j) = 0.5d0*( p(n_p,j) + p(n_p+1,j) )
        enddo
      endif

c--------------end interp lambda--------------------


      do j = 1, ny_n
c--------------copy u_e, v_e, rho_e to NS boundary----------------------
        re_bc_ew(j,1) = 0.5d0*
     &                   ( r_e(nx_start,ny_start+j-1) +
     &                     r_e(nx_start-1,ny_start+j-1) )
 
        vel1 = v_e(nx_start,   ny_start+j-1 )
        vel2 = v_e(nx_start,   ny_start+j   )
        vel3 = v_e(nx_start+1, ny_start+j-1 )
        vel4 = v_e(nx_start+1, ny_start+j   )

        v_bc_ew(j,1) = 0.25*(vel1+vel2+vel3+vel4)

        u_bc_ew(j,1) = u_e(nx_start, ny_start+j-1)

c---copy u_e,v_e for sending into add_ke routine-------
         ue_bc_ew(j,1) = u_bc_ew(j,1)
         ve_bc_ew(j,1) = v_bc_ew(j,1)
        
c-------------------end-copy--------------------------------------------- 

        if (tang .and. project) then
          v_bc_ew(j,1) = v_bc_ew(j,1)*re_bc_ew(j,1)
        endif


        if (project) then

          u_bc_ew(j,1) = u_proj(1+n_p,j+n_p)

          if (n_p .ge. 1) then
          v_bc_ew(j,1) = v_bc_ew(j,1) -
     &                   (p_y(n_p+2 +j-1) - p_y(n_p +j-1))/(2.d0*dy_n)
          endif


        endif

        if (tang .and. project) then
          v_bc_ew(j,1) = v_bc_ew(j,1)/density
        endif


         if (ke_flux .and. neum) then

           au_bc_ew(j,1) = (re_bc_ew(j,1)/density)*
     &                    ( u_bc_ew(j,1)**2 + v_bc_ew(j,1)**2 )*
     &                      u_bc_ew(j,1)

           av_bc_ew(j,1) = (re_bc_ew(j,1)/density)*
     &                    ( u_bc_ew(j,1)**2 + v_bc_ew(j,1)**2 )*
     &                      v_bc_ew(j,1)

           u_bc_ew(j,1) = au_bc_ew(j,1)*
     &                  ( au_bc_ew(j,1)**2 + av_bc_ew(j,1)**2)** 
     &                  (-1.d0/3.d0)

           v_bc_ew(j,1) = av_bc_ew(j,1)*
     &                  ( au_bc_ew(j,1)**2 + av_bc_ew(j,1)**2)**
     &                  (-1.d0/3.d0)

         endif ! ke_flux and neum

         if (sq .and. neum) then

          u_bc_ew(j,1) = u_bc_ew(j,1)*sqrt(re_bc_ew(j,1)/density)
          v_bc_ew(j,1) = v_bc_ew(j,1)*sqrt(re_bc_ew(j,1)/density)

         endif ! sq and neum

         if (mass .and. neum) then

          u_bc_ew(j,1) = u_bc_ew(j,1)*(re_bc_ew(j,1)/density)
          v_bc_ew(j,1) = v_bc_ew(j,1)*(re_bc_ew(j,1)/density)

         endif ! mass and neum
   
      enddo ! j=1,ny_n


! east

c--------------interpolate lambd vertic---------------

      if (n_p .ge. 1) then
        do j = 1,ny_p
           p_y(j) = 0.5d0*( p(n_p+nx_n,j) + p(n_p+nx_n+1,j) )
        enddo
      endif

c--------------end interp lambda--------------------


      do j = 1,ny_n

c--------------copy u_e, v_e, rho_e to NS boundary----------------------
        re_bc_ew(j,2) = 0.5d0*
     &                   ( r_e(nx_start+nx_n,ny_start+j-1) +
     &                     r_e(nx_start+nx_n-1,ny_start+j-1) )
 
        vel1 = v_e(nx_start+nx_n  , ny_start+j-1 )
        vel2 = v_e(nx_start+nx_n  , ny_start+j   )
        vel3 = v_e(nx_start+nx_n-1, ny_start+j-1 )
        vel4 = v_e(nx_start+nx_n-1, ny_start+j   )

        v_bc_ew(j,2) = 0.25*(vel1+vel2+vel3+vel4)

        u_bc_ew(j,2) = u_e(nx_start+nx_n,ny_start+j-1)

c---copy u_e,v_e for sending into add_ke routine-------
         ue_bc_ew(j,2) = u_bc_ew(j,2)
         ve_bc_ew(j,2) = v_bc_ew(j,2)

c--------------------end-copy---------------------------------------------

        if (tang .and. project) then
          v_bc_ew(j,2) = v_bc_ew(j,2)*re_bc_ew(j,2)
        endif


        if (project) then
          u_bc_ew(j,2) = u_proj(n_p+nx_n+1,j+n_p)

          if (n_p .ge. 1) then
          v_bc_ew(j,2) = v_bc_ew(j,2) -
     &                ( p_y(n_p+2 +j-1) - p_y(n_p +j-1) ) / (2.d0*dy_n)
          endif


        endif

 
        if (tang .and. project) then
          v_bc_ew(j,2) = v_bc_ew(j,2)/density
        endif



      enddo ! j = 1,ny_n

      return
      end



