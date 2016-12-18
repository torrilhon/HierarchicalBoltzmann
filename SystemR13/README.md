# Regularized 13-Moment-Equations Implementation

Equations and problem data for the steady, linear R13 system using equations for
density, velocity, temperature, stress and heat flux, as well as higher 
order moments m, and R (delta vanishes in steady linear case). 
 
### Usage:
Compile with `make` in the current directory after adjusting the 
directories in `Makefile???.def` in top directory. Use `make clean` 
if other examples have been compiled before. 

Run with specifying a grid file, e.g., 

    ./SystemR13.out ../grids/TwoCirclesLargeH0p2.msh

This solves the system with default parameters and generates (appends) the 
file `MasterLog.txt` in the output folder with informations, like parameters, 
mesh size and errors.

Additional comand line parameters are:

    -tau value        //= real-valued positive relaxation time
    -order value      //= 1 (linear) or 2 (quadratic elements)
    -output value     //= 0 or 1, writes data file in output directory
    -mfile filename   //= alternative logfile

The output file contains the result as simple list of vertex coordinates and 
data (rho, vx, vy, theta, stress, qx, qy), where stress is the von-Mises stress 
invariant. This file can be visualized with the mathematica script `Visualization.nb` 
in the evaluation folder.

More adjustments to the equation and boundary data are possible in the files 
`SystemR13.cpp` and `EquationsR13.cpp` and require a re-compile.

The batch file `BatchRunTwoCirclesLarge` will run a set of simulations on different 
grids and with different parameters for empirical error analysis. The resulting logfile 
can be evaluated with the mathematica script `ErrorEvaluation.nb` in the 
evaluation folder.
 
