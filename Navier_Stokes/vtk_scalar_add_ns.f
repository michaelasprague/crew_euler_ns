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

      subroutine vtk_scalar_add_ns(numnodes,temp,
     &                  filename,scalarname,ifinal)

      implicit double precision (a-h,o-z)
      integer numnodes
      double precision temp(numnodes)

      character filename*(*)
      character scalarname*(*)
      character command*80

      open (unit=20, file=filename//'.vtk', status='unknown')

      write(20,*) 'SCALARS ',scalarname,' float 1'
      write(20,*) 'LOOKUP_TABLE default'

      do i = 1,numnodes
        write(20,*) temp(i)
      enddo

      call flush(20)

      if (ifinal .eq. 1) then
        command = '../asa/asa '//filename//'.vtk > tmp.vtk'
        call system(command)
        command = '/bin/mv -f tmp.vtk '//filename//'.vtk'
        call system(command)
      endif

      return
      end
