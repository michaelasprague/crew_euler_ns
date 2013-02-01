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

      subroutine vtk_scalar_out_ns(nodes,elements,numel,numnodes,temp,
     &                  filename,scalarname)

      implicit double precision (a-h,o-z)
      integer numnodes,numel
      integer elements(4,numel)
      double precision nodes(2,numnodes)
      double precision temp(numnodes)

      character filename*(*)
      character scalarname*(*)
c     character command*80

      open (unit=20, file=filename//'.vtk', status='unknown')

      write(20,*)'# vtk DataFile Version 2.0'

      write(20,*) 'Really cool data'

      write(20,*) 'ASCII'

      write(20,*) 'DATASET UNSTRUCTURED_GRID'

      write(20,*) 'POINTS', numnodes, 'float'

! find abs max 
      dmax = 0.d0 
      do i = 1,numnodes
        if (abs(temp(i)) .gt. dmax) then
          dmax = temp(i) 
        endif
      enddo


      do i = 1,numnodes
        write(20,*) nodes(1,i)
        write(20,*) nodes(2,i)
        write(20,*) 0.
        !write(20,*) 0.1 * temp(i) / dmax
        write(20,*) ' '
      enddo

      write(20,*) 'CELLS', numel, 5*numel

      do i = 1,numel
        write(20,*) '4', (elements(l,i)-1,l = 1,4)
      enddo

      write(20,*) 'CELL_TYPES', numel
      write(20,*) (9, i = 1,numel)

      write(20,*) 'POINT_DATA', numnodes
      write(20,*) 'SCALARS ',scalarname,' float 1'
      write(20,*) 'LOOKUP_TABLE default'

      do i = 1,numnodes
        if (abs(temp(i)) .gt. 1d-9) then
          write(20,*) temp(i)
        else
          write(20,*) 0.
        endif
      enddo

      call flush(20)

      return
      end
