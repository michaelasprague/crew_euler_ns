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

      subroutine vtk_vector_add(nodes,elements,numel,numnodes,u,v,
     &                  filename,vectorname,ifinal)

      implicit double precision (a-h,o-z)
      integer norder,numnodes,numel
      integer elements(4,numel)
      double precision nodes(2,numnodes)
      double precision u(numnodes)
      double precision v(numnodes)

      character filename*(*)
      character vectorname*(*)
      character command*80

      open (unit=20, file=filename//'.vtk', status='unknown')

      write(20,*) 'VECTORS ',vectorname,' float'

      do i = 1,numnodes
        write(20,*) u(i),v(i), 0.
      enddo

      call flush(20)

      close(unit=20)

      if (ifinal .eq. 1) then
        command = '../asa/asa '//filename//'.vtk > tmp.vtk'
        call system(command)
        command = '/bin/mv -f tmp.vtk '//filename//'.vtk'
        call system(command)
      endif

      return
      end
