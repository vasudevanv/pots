# Non-Standard Interaction Potentials 

This repository contains scripts and codes which enable the use of some
non-standard interaction potentials with molecular simulation software, 
mainly gromacs and lammps. 

## WCA scaling - 12-6 potential

The Weeks-Chandler-Anderson (WCA) model potential is a purely repulsive
potential used to model excluded volume interactions. The potential is 
constructed by truncating the Lennard-Jones potential at r = 2^{1/6} \sigma
and shifting it by \epsilon ensuring that the potential decays smoothly to 0
at r = 2^{1/6} \sigma. This scheme can be used more generally to scale the 
attractive portion of the interaction potential without much modification to 
the repulsive portion. A scaling value of 1.0 results in the normal 
Lennard-Jones potential. The scaled potential is given by:

           4*\epsilon*[ (sig/r)**12 - (sig/r)**6 + (1 - \lambda)/4 ] , r< 2**(1/6) \sigma 
    u(r) =
           4*\epsilon*\lambda*[ (sig/r)**12 - (sig/r)**6 ] , r >= 2**(1/6) \sigma 
           
### GROMACS
The gromacs tabulated potential option provides an easy way to use the WCA 
scaled potential. The table-gmx.f90/table-gmx.py codes produce interaction 
tables for use with gromacs. 

### LAMMPS
For use with lammps, we have implemented new pair_style's based on pair_style 
lj/cut/coul/long and pair_style lj/cut. The potential can be used with the 
following syntax:

* style name:    lj/cut/coul/long/wca 
* options:       pair_style lj/cut/coul/long/wca <cut_lj> <lambda> <cut_coulomb>
* coefficients:  pair_coeff <I> <J> <epsilon> <sigma> <lambda> <cut_lj>

* style name:    lj/cut/wca 
* options:       pair_style lj/cut/wca <cut_lj> <lambda>
* coefficients:  pair_coeff <I> <J> <epsilon> <sigma> <lambda> <cut_lj_one>

