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

      subroutine Hbar_vec(w, u, ew_rc, ns_rc,
     &                    nx, ny, dx, dy, dt, viscosity)


      implicit double precision (a-h, o-z)

      double precision w(nx, ny)
      double precision u(nx, ny)

C these are u_n and v_n interpolated to the boundaries needed for H^n
      double precision ew_rc(nx+1, ny)
      double precision ns_rc(nx,   ny+1)

! the divide-by-two comes from averaging

C inside
c!$omp parallel private(i,j) 
c!$omp& shared(w,u,ew_rc,ns_rc,dy,dx,viscosity,nx,ny)
c!$omp do
      do j = 2, ny-1

        do i = 2,nx-1

          w(i,j) =  
     &              
     &           u(i-1,j)* ( ew_rc(i,j)  * dy/2.d0 + viscosity*dy/dx )
     &
     &         + u(i+1,j)* (-ew_rc(i+1,j)* dy/2.d0 + viscosity*dy/dx )
     &
     &         + u(i,j-1)* ( ns_rc(i,j)  * dx/2.d0 + viscosity*dx/dy )
     &
     &         + u(i,j+1)* (-ns_rc(i,j+1)* dx/2.d0 + viscosity*dx/dy )

!mas checks 9/9/2011           
 
        enddo

      enddo
c!$omp end parallel

C West


      do j = 2, ny-1

        i = 1
        
         w(i,j) =   
     &                       
     &          0.d0 
     &
     &         + u(i+1,j)* (-ew_rc(i+1,j)* dy/2.d0 + viscosity*dy/dx )
     &
     &         + u(i,j-1)* ( ns_rc(i,j)  * dx/2.d0 + viscosity*dx/dy )
     &
     &         + u(i,j+1)* (-ns_rc(i,j+1)* dx/2.d0 + viscosity*dx/dy )

C East

        i = nx

         w(i,j) =  
     &
     &          u(i-1,j)* ( ew_rc(i,j)  * dy/2.d0 + viscosity*dy/dx )
     &
     &         + 0.d0 
     &
     &         + u(i,j-1)* ( ns_rc(i,j)  * dx/2.d0 + viscosity*dx/dy )
     &
     &         + u(i,j+1)* (-ns_rc(i,j+1)* dx/2.d0 + viscosity*dx/dy )

      enddo


C South

      do i = 2,nx-1

        j = 1

         w(i,j) =  
     &
     &           u(i-1,j)* ( ew_rc(i,j)  * dy/2.d0 + viscosity*dy/dx )
     &
     &         + u(i+1,j)* (-ew_rc(i+1,j)* dy/2.d0 + viscosity*dy/dx )
     &
     &         + 0.d0 
     &
     &         + u(i,j+1)* (-ns_rc(i,j+1)* dx/2.d0 + viscosity*dx/dy )

C North

        j = ny

          w(i,j) =  
     &
     &           u(i-1,j)* ( ew_rc(i,j)  * dy/2.d0 + viscosity*dy/dx )
     &
     &         + u(i+1,j)* (-ew_rc(i+1,j)* dy/2.d0 + viscosity*dy/dx )
     &
     &         + u(i,j-1)* ( ns_rc(i,j)  * dx/2.d0 + viscosity*dx/dy )
     &
     &         + 0.d0 

      enddo

C corners:

C SW

      j = 1
      i = 1

        w(i,j) =  
     &
     &          0.d0 
     &
     &         + u(i+1,j)* (-ew_rc(i+1,j)* dy/2.d0 + viscosity*dy/dx )
     &
     &         + 0.d0 
     &
     &         + u(i,j+1)* (-ns_rc(i,j+1)* dx/2.d0 + viscosity*dx/dy )


C SE

      j = 1
      i = nx

        w(i,j) =  
     &
     &          u(i-1,j)* ( ew_rc(i,j)  * dy/2.d0 + viscosity*dy/dx )
     &
     &         + 0.d0 
     &
     &         + 0.d0  
     &
     &         + u(i,j+1)* (-ns_rc(i,j+1)* dx/2.d0 + viscosity*dx/dy )


C NE

      j = ny
      i = nx

         w(i,j) =  
     &
     &          u(i-1,j)* ( ew_rc(i,j)  * dy/2.d0 + viscosity*dy/dx )
     &
     &         + 0.d0 
     &
     &         + u(i,j-1)* ( ns_rc(i,j)  * dx/2.d0 + viscosity*dx/dy )
     &
     &         + 0.d0 

C NW

      j = ny
      i = 1

         w(i,j) =  
     &
     &          0.d0 
     &
     &         + u(i+1,j)* (-ew_rc(i+1,j)* dy/2.d0 + viscosity*dy/dx )
     &
     &         + u(i,j-1)* ( ns_rc(i,j)  * dx/2.d0 + viscosity*dx/dy )
     &
     &         + 0.d0 


C DONE!


      return
      end
 


       
                             
 

       


