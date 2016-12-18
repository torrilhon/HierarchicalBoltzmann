# System A Implementation

Equations, problem data and makefile for model system (A). 
 
### Usage:
Compile with `make` in the current directory after adjusting the 
directories in `Makefile???.def` in top directory. Use `make clean` 
if other examples have been compiled before. 

Run with specifying a grid file, e.g., 

    ./SystemA.out ../grids/TwoCirclesLargeH0p2.msh

This solves the system with default parameters and generates (appends) the 
file `MasterLog.txt` in the output folder with informations, like parameters, 
mesh size and errors.

Additional comand line parameters are:

    -tau value        //= real-valued positive relaxation time
    -order value      //= 1 (linear) or 2 (quadratic elements)
    -kappa value      //= real-valued parameter in boundary conditions
    -output value     //= 0 or 1, writes data file in output directory
    -mfile filename   //= alternative logfile

For the role of `kappa` please consult the description in [MT2016].
The output file contains the result as simple list of vertex coordinates and 
data (theta, qx, qy). This file can be visualized with the mathematica script
`Visualization.nb` in the evaluation folder.

More adjustments to the equation and boundary data are possible in the files 
`SystemA.cpp` and `EquationsA.cpp` and require a re-compile.

The batch file `BatchRun` will run a set of simulations on different grids and 
with different parameters for empirical error analysis. The resulting logfile 
can be evaluated with the mathematica script `ErrorEvaluation.nb` in the 
evaluation folder.
