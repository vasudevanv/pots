# Tabulated potentials

This repository contains codes which generate interaction potential tables 
for use with molecular simulation software, mainly gromacs.

## WCA scaling - 12-6 potential

The Weeks-Chandler-Anderson (WCA) model potential is a purely repulsive
potential used to model excluded volume interactions. The potential is 
constructed by truncating the Lennard-Jones potential at r = 2^{1/6} \sigma
and shifting it by \epsilon ensuring that the potential decays smoothly to 0
at r = 2^{1/6} \sigma.

           4*\epsilon*[ (sig/r)**12 - (sig/r)**6 + (1 - \lambda)/4 ] , r< 2**(1/6) \sigma 
    u(r) =
           4*\epsilon*\lambda*[ (sig/r)**12 - (sig/r)**6 ] , r >= 2**(1/6) \sigma 
           
This scheme can be used more generally to scale the attractive interactions 
without much modification to the repulsive portion. A scaling value of 1.0 
results in the normal Lennard-Jones potential.

### GROMACS
The table-gmx.f90/table-gmx.py codes produce interaction tables for interactions scaled 
using the WCA scheme for use with gromacs. 

### LAMMPS
A LAMMPS pair potential has also been implemented based on pair_style lj/cut/coul/long
and pair_style lj/cut. The potential can be used with the following syntax:

* style name:    lj/cut/coul/long/wca 
* options:       pair_style lj/cut/coul/long/wca <cut_lj> <lambda> <cut_coulomb>
* coefficients:  pair_coeff <I> <J> <epsilon> <sigma> <lambda> <cut_lj>

* style name:    lj/cut/wca 
* options:       pair_style lj/cut/wca <cut_lj> <lambda>
* coefficients:  pair_coeff <I> <J> <epsilon> <sigma> <lambda> <cut_lj_one>

