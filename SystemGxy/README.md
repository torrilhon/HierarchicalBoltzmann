# Hierarchical Boltzmann Framework

Equations and problem data for hierarchical Boltzmann simulation based on 
Hermite-discretizations (steady, linear, Maxwell molecules).
 
### Usage:
Compile with `make` in the current directory after adjusting the 
directories in `Makefile???.def` in top directory. Use `make clean` 
if other examples have been compiled before. 

Run with specifying a grid file, e.g., 

    ./SystemGxy.out ../grids/complexH0p04.msh

This solves the system with default parameters on the given mesh and 
generates (appends) the file `MasterLog.txt` in the output folder with 
informations, like parameters, mesh size, etc.

Additional command line parameters are:

    -moments value    //= integer value for number of moments used
    -tau value        //= real-valued positive relaxation time
    -order value      //= 1 (linear) or 2 (quadratic elements)
    -output value     //= 0 or 1, writes data file in output directory
    -mfile filename   //= alternative logfile

The output file contains the result as simple list of vertex coordinates and 
data (rho, vx, vy, theta, stress, qx, qy), where stress is the von-Mises stress 
invariant. This file can be visualized with the mathematica script `Visualization.nb` 
in the evaluation folder.

More adjustments to the equation and boundary data are possible in the files 
`SystemGxy.cpp` and `EquationsGxy.cpp` and require a re-compile.

The batch file `BatchRunComplexFlow` will run a set of simulations with 
different models for model error estimation.

