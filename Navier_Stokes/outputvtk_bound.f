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

      subroutine outputvtk_bound(u_bc_ew,u_bc_ns,v_bc_ew,v_bc_ns,
     &           xmin,ymin,xmax,ymax, umean, vmean,
     &           dx,dy,nx,ny,filename)

      implicit double precision (a-h,o-z)

      double precision u_bc_ew(ny,2), u_bc_ns(nx,2)
      double precision v_bc_ew(ny,2), v_bc_ns(nx,2)

      integer nelem(2,2*(nx-1)+2*(ny-1))
 
      integer ihat_array(4)

      double precision xyz(2,2*nx + 2*ny)
      double precision u_vtk(2*nx + 2*ny)
      double precision v_vtk(2*nx + 2*ny)

      character filename*(*)

      character command*80

      nmax = (nx)*(ny) 
      nemax = (nx-1)*(ny-1)

      y = ymin 
      x = xmin + 0.5*dx
      ihat = 0
      do i = 1, nx

        ihat =  ihat + 1 

        xyz(1,ihat) = x
        xyz(2,ihat) = y 

        u_vtk(ihat) = u_bc_ns(i,1)
        v_vtk(ihat) = v_bc_ns(i,1)

! subtract mean
        u_vtk(ihat) = u_vtk(ihat) - umean
        v_vtk(ihat) = v_vtk(ihat) - vmean

        x = x + dx
      enddo

      y = ymin + 0.5*dy
      x = xmax
      do i = 1, ny

        ihat =  ihat + 1 

        xyz(1,ihat) = x
        xyz(2,ihat) = y 

        u_vtk(ihat) = u_bc_ew(i,2)
        v_vtk(ihat) = v_bc_ew(i,2)

! subtract mean
        u_vtk(ihat) = u_vtk(ihat) - umean
        v_vtk(ihat) = v_vtk(ihat) - vmean

        y = y + dy
      enddo

      y = ymax
      x = xmin + 0.5*dx
      do i = 1, nx

        ihat =  ihat + 1 

        xyz(1,ihat) = x
        xyz(2,ihat) = y 

        u_vtk(ihat) = u_bc_ns(i,2)
        v_vtk(ihat) = v_bc_ns(i,2)

! subtract mean
        u_vtk(ihat) = u_vtk(ihat) - umean
        v_vtk(ihat) = v_vtk(ihat) - vmean

        x = x + dx
      enddo

      y = ymin + 0.5*dy
      x = xmin
      do i = 1, ny
        ihat =  ihat + 1 

        xyz(1,ihat) = x
        xyz(2,ihat) = y 

        u_vtk(ihat) = u_bc_ew(i,1)
        v_vtk(ihat) = v_bc_ew(i,1)

! subtract mean
        u_vtk(ihat) = u_vtk(ihat) - umean
        v_vtk(ihat) = v_vtk(ihat) - vmean

        y = y + dy
      enddo


      ihat_array(1) = 0
      ihat_array(2) = nx
      ihat_array(3) = nx + ny
      ihat_array(4) = nx + ny + nx

      ielem = 0
      do j = 1, 4
        ihat = ihat_array(j)
        do i = 1, nx -1
          ihat = ihat + 1
          ielem = ielem + 1 
          nelem(1,ielem) = ihat
          nelem(2,ielem) = ihat + 1
        enddo
      enddo

      numel = ielem

      numnodes = 2*nx + 2*ny


      open (unit=20, file=filename//'.vtk', status='unknown')


      write(20,*)'# vtk DataFile Version 2.0'

      write(20,*) 'Really cool data'

      write(20,*) 'ASCII'

      write(20,*) 'DATASET UNSTRUCTURED_GRID'

      write(20,*) 'POINTS', numnodes, 'float'

      do i = 1,numnodes
        write(20,991) xyz(1,i)
        write(20,991) xyz(2,i)
        write(20,991) 0.
        !write(20,*) 0.1 * temp(i) / dmax
        write(20,*) ' '
      enddo


      write(20,*) 'CELLS', numel, 3*numel

      do i = 1,numel
        write(20,*) '2', (nelem(l,i)-1,l = 1,2)
      enddo

      write(20,*) 'CELL_TYPES', numel
      write(20,*) (3, i = 1,numel)

      write(20,*) 'POINT_DATA', numnodes
      write(20,*) 'SCALARS ','zero',' float 1'
      write(20,*) 'LOOKUP_TABLE default'

      do i = 1,numnodes
         write(20,991) 0.
      enddo

      write(20,*) 'VECTORS ','boundvectors',' float'

      do i = 1,numnodes
        uout = u_vtk(i)
        vout = v_vtk(i)
        if (abs(uout) .lt. 1d-10) uout = 0.d0
        if (abs(vout) .lt. 1d-10) vout = 0.d0
        write(20,993) uout,vout, 0.
      enddo


      call flush(20)

      close(unit=20)

      ifinal = 1

      if (ifinal .eq. 1) then
        !command = '/usr/bin/asa '//filename//'.vtk > tmp.vtk'
        command = '../asa/asa '//filename//'.vtk > tmp.vtk'
        call system(command)
        command = '/bin/mv -f tmp.vtk '//filename//'.vtk'
        call system(command)
      endif

991   format(E16.8)
993   format(3E16.8)

      return
      end

