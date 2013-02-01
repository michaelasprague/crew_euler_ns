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

11,12c11,12
<       parameter(nx_e_quarter = 50) !quarter size of Euler (must be even) 
<       parameter(ny_e_quarter = 50) !quarter size of Euler (must be even) 
---
>       parameter(nx_e_quarter = 80) !quarter size of Euler (must be even) 
>       parameter(ny_e_quarter = 80) !quarter size of Euler (must be even) 
100c100
<       double precision energy2d_n0  (nx_n, ny_n) !KE density  
---
>       double precision energy2d0_n  (nx_n, ny_n) !KE density  
227c227
<       project = .true.
---
>       project = .false.
237c237
<       neum = .false. ! Neumann bcs for NS domain (North and East)
---
>       neum = .true. ! Neumann bcs for NS domain (North and East)
253c253
<       ke_add  = .true. ! option for  proj
---
>       ke_add  = .false. ! option for  proj
771,781c771,781
<         b_flux = 0.d0
<         do i = 1,ny_n
<            b_flux = b_flux - u_bc_ew(i,1)*dy_n + u_bc_ew(i,2)*dy_n
<         enddo
<         do i = 1,nx_n
<            b_flux = b_flux - v_bc_ns(i,1)*dx_n + v_bc_ns(i,2)*dx_n
<         enddo
<         if (abs(b_flux) .gt. 10.*tol_cg_proj) then
<            write(*,*) 'b_flux exceeds limit:', b_flux
<            stop ' b_flux exceeds limit'
<         endif
---
> c        b_flux = 0.d0
> c        do i = 1,ny_n
> c           b_flux = b_flux - u_bc_ew(i,1)*dy_n + u_bc_ew(i,2)*dy_n
> c        enddo
> c        do i = 1,nx_n
> c           b_flux = b_flux - v_bc_ns(i,1)*dx_n + v_bc_ns(i,2)*dx_n
> c        enddo
> c        if (abs(b_flux) .gt. 10.*tol_cg_proj) then
> c           write(*,*) 'b_flux exceeds limit:', b_flux
> c           stop ' b_flux exceeds limit'
> c        endif
899c899
<           call calc_energy_n(energy_n0, energy2d_n0, u0_n, v0_n, r0_n,
---
>           call calc_energy_n(energy0_n, energy2d0_n, u0_n, v0_n, r0_n,
900a901
> 
903c904
<           call error_ns_ke(energy2d_n, energy2d_n0, err_energy2d_n,
---
>           call error_ns_ke(energy2d_n, energy2d0_n, err_energy2d_n,
908c909
<           write(*,*)'energy', i_vtk_out, t, energy_n, energy_n0, alpha
---
>           write(*,*)'energy', i_vtk_out, t, energy_n, energy0_n, alpha
911c912
<           write(54,*) i_vtk_out, t, energy_n, energy_n0, alpha
---
>           write(54,*) i_vtk_out, t, energy_n, energy0_n, alpha
982a984,987
> c ENERGY0
> 
>             call matlabvtk(energy2d0_n, nx_n, ny_n,
>      &              'energy2d0_n'//filename(i_vtk_out))
988a994,998
> c NS0
>             call matlabvtk(u0_n,nx_n,ny_n,
>      &              'u0_vel_n'//filename(i_vtk_out))
>             call matlabvtk(v0_n,nx_n,ny_n,
>      &               'v0_vel_n'//filename(i_vtk_out))
