#!/usr/bin/env python
#
#  Python  code to calculate VdW lookup tables for Gromacs 4.0.X
# 
#  Copyright 2012, Vasudevan Venkateshwaran
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
#  table.py
# 
#     Usage:
#     python table.py GROUPA SIGA EPSA GROUPB SIGB EPSB SCALE (VERSION)
#     python table.py SOL 0.30 0.60 ATM 0.44 0.85 0.50 (4)
#     
#     Version is by default 4. If you want tables for gromacs 3.X.X use
#     version = 3. The default version = 4 produces tables which can be
#     used with gromacs 4.X.X

import math
import string
import argparse

def wcatables(groupa,siga,epsa,groupb,sigb,epsb,scale,version=4):
    
    # Table Parameters
    xmax = 10.0              # in nm
    dx = 0.002               # in nm
    bmax = int(xmax/dx)+1
  
    # Create the table
    fname="table_" + groupa.strip() + "_" + groupb.strip() + ".xvg"
    eps =  math.sqrt(epsa*epsb)
    sigma = (siga + sigb)/2.0
    csix = 4*eps*sigma**6
    ctwel = 4*eps*sigma**12
    f = open(fname,'w')

    # Set output format
    fmt = "%14.5e" * 7 + "\n"

    # Generate the table
    for i in range(bmax+1):
        x = i * dx
        if (x == 0) :  
            f.write( fmt % (0.0,0.0,0.0,0.0,0.0,0.0,0.0) )
        else:
            if(x <= 2**(1.0/6.0)*sigma) :
                attx = -1.0/x**6.
                attxd = 6.0/x**7
                attxdd = -42.0/x**8
                repx = 1.0/x**12 + (1-scale)/(4.0*sigma**12)
                repxd = -12.0/x**13
                repxdd = 156.0/x**14
                if (version == 4) :
                    f.write( fmt %
                             (x,1.0/x,1.0/x**2,attx,-attxd,repx,-repxd) )
                elif (version == 3) :
                    f.write( fmt % 
                             (x,1.0/x,1.0/x**2,attx,attxdd,repx,repxdd) )
            else:
                attx = -scale/x**6.
                attxd = 6.0*scale/x**7
                attxdd = -42.0*scale/x**8
                repx = scale*1.0/x**12
                repxd = -12.0*scale*1.0/x**13
                repxdd = scale*156.0/x**14
                if (version == 4) :
                    f.write( fmt % 
                             (x,1.0/x,1.0/x**2,attx,-attxd,repx,-repxd) )
                elif (version == 3) :
                    f.write( fmt % 
                             (x,1.0/x,1.0/x**2,attx,attxdd,repx,repxdd) )
                    
        
    f.close()
  

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # Group 1 parameters
    parser.add_argument('groupa', action='store',
                        dest='groupa',
                        help='Name of group 1')
    parser.add_argument('siga', action='store', 
                        dest='boxx', type=float, 
                        help='LJ Sigma for group 1')
    parser.add_argument('epsa', action='store', 
                        dest='epsa', type=float, 
                        help='LJ Epsilon for group 1')
    # Group 2 parameters
    parser.add_argument('groupb', action='store',
                        dest='groupb',
                        help='Name of group 2')
    parser.add_argument('sigb', action='store', 
                        dest='sigb', type=float, 
                        help='LJ Sigma for group 2')
    parser.add_argument('epsb', action='store', 
                        dest='epsb', type=float, 
                        help='LJ Epsilon for group 2')
    # Scaling
    parser.add_argument('scale', action='store', 
                        dest='scale', type=float, 
                        help='WCA scaling parameter')
    # Gromacs version
    parser.add_argument('version', action='store', 
                        dest='version', type=int, default=4,
                        help='Gromacs version (3/4). Default is 4')
    # Parse arguments
    args = parser.parse_args()
    # Call function
    wcatables(args.groupa, args.siga, args.epsa,
              args.groupb, args.sigb, args.epsb,
              scale, version)


  
