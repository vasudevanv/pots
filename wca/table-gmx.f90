! Fortran 90 code to calculate VdW lookup tables for Gromacs
!
! Copyright 2012, Vasudevan Venkateshwaran
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! table.f90
!
!    Compile : gfortran table.f90 -o gentable
!
!    Usage:
!    ./gentable GROUPA SIGA EPSA GROUPB SIGB EPSB SCALE (VERSION)
!    ./gentable SOL 0.30 0.60 ATM 0.44 0.85 0.50 (4)
!    
!    Version is by default 4. If you want tables for gromacs 3.X.X use
!    version = 3. The default version = 4 produces tables which can be
!    used with gromacs 4.X.X

program wcatables
  implicit none
  integer         i,j,k,bmax
  integer         version
  real*8          siga,epsa,sigb,epsb
  real*8          x,xmax,dx,scale,sqscale,sigma,eps,csix,ctwel
  real*8          attx,repx,attxd,repxd,attxdd,repxdd
  character*14    csiga,cepsa,csigb,cepsb,cscale,cversion
  character*5     groupa,groupb
  character*80    fname
  
  ! Table Parameters
  xmax = 10                ! in nm
  dx = 0.002               ! in nm
  bmax = int(xmax/dx)+1
  
  ! Read input arguments
  if((iargc().eq.7).or.(iargc().eq.8)) then
     call getarg(1,groupa)
     call getarg(2,csiga)
     call getarg(3,cepsa)
     call getarg(4,groupb)
     call getarg(5,csigb)
     call getarg(6,cepsb)
     call getarg(7,cscale)
     read(csiga,*) siga
     read(cepsa,*) epsa
     read(csigb,*) sigb
     read(cepsb,*) epsb
     read(cscale,*) scale
     if (iargc().eq.8) then
        call getarg(8,cversion)
        read(cversion,*) version
     else
        version = 4
     endif
     write(*,*) ' Scaling ',groupa,'-',groupb,' attractions by : ', &
          scale
  else
     write(*,*) "Usage: group1 sig1 eps1 group2 sig2 eps2 scale (version)"
     write(*,*) "Note : version is by default 4"
  endif
  
  ! Create the table
  fname="table_"//trim(groupa)//"_"//trim(groupb)//".xvg"
  eps =  sqrt(epsa*epsb)
  sigma = (siga + sigb)/2.0
  csix = 4*eps*sigma**6
  ctwel = 4*eps*sigma**12
  open(unit = 24, file = fname)
  do i = 0, bmax
     x = i * dx
     if(x.eq.0)then 
        write(24,'(7(1pe14.5))')0.0,0.0,0.0,0.0,0.0,0.0,0.0
     else
        if(x.le.2**(1.0/6.0)*sigma) then
           attx = -1.0/x**6.
           attxd = 6.0/x**7
           attxdd = -42.0/x**8
           repx = 1.0/x**12 + (1-scale)/(4.0*sigma**12)
           repxd = -12.0/x**13
           repxdd = 156.0/x**14
           if (version.eq.4) then
              write(24,'(7(1pe14.5))') &
                   x,1.0/x,1.0/x**2,attx,-attxd,repx,-repxd
           else if (version.eq.3) then
              write(24,'(7(1pe14.5))') &
                   x,1.0/x,1.0/x**2,attx,attxdd,repx,repxdd
           endif
        else
           attx = -scale/x**6.
           attxd = 6.0*scale/x**7
           attxdd = -42.0*scale/x**8
           repx = scale*1./x**12
           repxd = -12.0*scale*1./x**13
           repxdd = scale*156.0/x**14
           if (version.eq.4) then
              write(24,'(7(1pe14.5))') &
                   x,1.0/x,1.0/x**2,attx,-attxd,repx,-repxd
           else if (version.eq.3) then
              write(24,'(7(1pe14.5))') &
                   x,1.0/x,1.0/x**2,attx,attxdd,repx,repxdd
           endif
        endif
     endif
  enddo
  close(24)
  
  stop
end program wcatables

