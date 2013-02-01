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

      subroutine A_vec(w, u, ew_rc, ns_rc,
     &                 nx, ny, dx, dy, dt, viscosity, neum, conv)

c diagonal part of H operator: A times u


      implicit double precision (a-h, o-z)

      double precision w(nx,ny), u(nx,ny)

C these are u_n and v_n interpolated to the boundaries needed for H^n
      double precision ew_rc(nx+1, ny)
      double precision ns_rc(nx , ny+1)

      logical neum
      logical conv 

C inside

C!$omp parallel private(i,j) 
C!$omp& shared(w,u,ew_rc,ns_rc,viscosity,dy,dx,nx,ny)
C!$omp do
      do j = 2, ny-1

        do i = 2,nx-1

          w(i,j) = u(i,j) 
     &            * ( ew_rc(i,j)*dy/2.d0 - ew_rc(i+1,j)*dy/2.d0
     &               +ns_rc(i,j)*dx/2.d0 - ns_rc(i,j+1)*dx/2.d0
     &               -viscosity *(dy/dx + dy/dx + dx/dy + dx/dy) )
            
        enddo

      enddo
C!$omp end parallel

C West
      do j = 2, ny-1

        i = 1
        
         w(i,j) = u(i,j)  
     &               *(        0.d0          - ew_rc(i+1,j)*dy/2.d0
     &               +ns_rc(i,j)*dx/2.d0 - ns_rc(i,j+1)*dx/2.d0
     &               -viscosity *(dy/dx +2.d0* dy/dx + dx/dy + dx/dy) )

! Why is there the 2 on the one dy/dx?  Has to do with evaluating the
! derivative at the boundary there.  Using centered value and boundary
! value; delta x is just dx/2 instead of dx  -- make sure factor of 2 is
! also in bn and A_vec 

C East

      i = nx

      if (neum) then

         w(i,j) = u(i,j)
     &             *( ew_rc(i,j)*dy/2.d0 - ew_rc(i+1,j)*dy 
     &               +ns_rc(i,j)*dx/2.d0 - ns_rc(i,j+1)*dx/2.d0
     &               -viscosity *(0.d0 + dy/dx + dx/dy + dx/dy))

      else if (conv) then

        term_e = (2.d0*dt*dy*ew_rc(i+1,j))/(dx+2.d0*dt*ew_rc(i+1,j))
        w(i,j) = u(i,j)
     &            *( ew_rc(i,j)*dy/2.d0 - ew_rc(i+1,j)* term_e 
     &               +ns_rc(i,j)*dx/2.d0 - ns_rc(i,j+1)*dx/2.d0
     &               -viscosity *(2.d0*dy/(dx+ 2.d0*dt*ew_rc(i+1,j)) 
     &                            + dy/dx + dx/dy + dx/dy))
 
      else

         w(i,j) = u(i,j) 
     &             *( ew_rc(i,j)*dy/2.d0 -         0.d0 
     &               +ns_rc(i,j)*dx/2.d0 - ns_rc(i,j+1)*dx/2.d0
     &               -viscosity *(2.d0*dy/dx + dy/dx + dx/dy + dx/dy))

      endif

      enddo


C South

      do i = 2,nx-1

        j = 1

         w(i,j) = u(i,j) 
     &            *( ew_rc(i,j)*dy/2.d0 - ew_rc(i+1,j)*dy/2.d0
     &               +      0.d0           - ns_rc(i,j+1)*dx/2.d0
     &           -viscosity *(dy/dx + dy/dx + dx/dy +2.d0*dx/dy) )

C North 

        j = ny

        if (neum) then

          w(i,j) = u(i,j)
     &                  *( ew_rc(i,j)*dy/2.d0 - ew_rc(i+1,j)*dy/2.d0
     &                    +ns_rc(i,j)*dx/2.d0 - ns_rc(i,j+1)*dx 
     &                -viscosity *(dy/dx + dy/dx +0.d0 + dx/dy) )
        
        else if (conv) then

        term_n = (2.d0*dt*dx*ns_rc(i,j+1))/(dy+2.d0*dt*ns_rc(i,j+1))
         w(i,j) = u(i,j)
     &             *( ew_rc(i,j)*dy/2.d0 - ew_rc(i+1,j)*dy/2.d0 
     &               +ns_rc(i,j)*dx/2.d0 - ns_rc(i,j+1)*term_n
     &               -viscosity *(2.d0*dx/(dy+ 2.d0*dt*ew_rc(i,j+1))
     &                            + dy/dx + dy/dx + dx/dy))

        else

          w(i,j) = u(i,j) 
     &                  *( ew_rc(i,j)*dy/2.d0 - ew_rc(i+1,j)*dy/2.d0
     &                    +ns_rc(i,j)*dx/2.d0 -        0.d0
     &                -viscosity *(dy/dx + dy/dx +2.d0*dx/dy + dx/dy) )
        endif

      enddo

C corners:

C SW

      j = 1
      i = 1

        w(i,j) = u(i,j) 
     &               *(       0.d0           - ew_rc(i+1,j)*dy/2.d0
     &               +       0.d0          - ns_rc(i,j+1)*dx/2.d0
     &        -viscosity *(dy/dx + 2.d0*dy/dx + dx/dy + 2.d0*dx/dy) )


C SE

      j = 1
      i = nx

! mas checks 2012-02-10

      if (neum) then

         w(i,j) = u(i,j)
     &        * ( ew_rc(i,j)*dy/2.d0 - ew_rc(i+1,j)*dy
     &        +        0.d0         - ns_rc(i,j+1)*dx/2.d0
     &        -viscosity *(0.d0 + dy/dx + dx/dy + 2.d0*dx/dy) )

      else if (conv) then

        term_e = (2.d0*dt*dy*ew_rc(i+1,j))/(dx+2.d0*dt*ew_rc(i+1,j))
         w(i,j) = u(i,j)
     &             *( ew_rc(i,j)*dy/2.d0 - ew_rc(i+1,j) * term_e
     &               +   0.d0            - ns_rc(i,j+1)*dx/2.d0
     &               -viscosity *(2.d0*dy/(dx+ 2.d0*dt*ew_rc(i+1,j))
     &                            + 3.d0*dx/dy + dy/dx))

      else

         w(i,j) = u(i,j) 
     &        * ( ew_rc(i,j)*dy/2.d0 -         0.d0
     &        +        0.d0         - ns_rc(i,j+1)*dx/2.d0
     &        -viscosity *(2.d0*dy/dx + dy/dx + dx/dy + 2.d0*dx/dy) )

      endif

C NE

      j = ny
      i = nx

      if (neum) then
         
         w(i,j) = u(i,j)
     &                  *( ew_rc(i,j)*dy/2.d0 - ew_rc(i+1,j)*dy 
     &                    +ns_rc(i,j)*dx/2.d0 - ns_rc(i,j+1)*dx 
     &       -viscosity *(0.d0 + dy/dx + 0.d0 + dx/dy) )

      else if (conv) then

        term_e = (2.d0*dt*dy*ew_rc(i+1,j))/(dx+2.d0*dt*ew_rc(i+1,j)) 
        term_n = (2.d0*dt*dx*ns_rc(i,j+1))/(dy+2.d0*dt*ns_rc(i,j+1)) 
         w(i,j) = u(i,j)
     &             *( ew_rc(i,j)*dy/2.d0 - ew_rc(i+1,j)*term_e
     &               +ew_rc(i,j)*dx/2.d0 - ns_rc(i,j+1)*term_n
     &               -viscosity *(2.d0*dy/(dx+ 2.d0*dt*ew_rc(i+1,j))+
     &                            2.d0*dx/(dy+ 2.d0*dt*ns_rc(i,j+1))
     &                            + dx/dy + dy/dx))

      else

         w(i,j) = u(i,j) 
     &                  *( ew_rc(i,j)*dy/2.d0 -         0.d0
     &                    +ns_rc(i,j)*dx/2.d0 -         0.d0
     &       -viscosity *(2.d0*dy/dx + dy/dx + 2.d0*dx/dy + dx/dy) )
 
      endif

C NW

      j = ny
      i = 1

      if (neum) then

         w(i,j) = u(i,j)
     &              *(         0.d0         - ew_rc(i+1,j)*dy/2.d0
     &                  +ns_rc(i,j)*dx/2.d0 - ns_rc(i,j+1)*dx 
     &        -viscosity *(dy/dx + 2.d0*dy/dx + 0.d0 + dx/dy) )

       else if (conv) then

        term_n = (2.d0*dt*dx*ns_rc(i,j+1))/(dy+2.d0*dt*ns_rc(i,j+1))
         w(i,j) = u(i,j)
     &             *( 0.d0               - ew_rc(i+1,j)*term_e
     &               +ew_rc(i,j)*dx/2.d0 - ns_rc(i,j+1)*dx/2.d0
     &               -viscosity *(2.d0*dx/(dy+ 2.d0*dt*ns_rc(i,j+1))
     &                            + dx/dy + 3.d0*dy/dx))

      else
 
 
         w(i,j) = u(i,j) 
     &              *(         0.d0         - ew_rc(i+1,j)*dy/2.d0
     &                  +ns_rc(i,j)*dx/2.d0 -          0.d0
     &        -viscosity *(dy/dx + 2.d0*dy/dx + 2.d0*dx/dy + dx/dy) )

      endif

C DONE!


      return
      end
