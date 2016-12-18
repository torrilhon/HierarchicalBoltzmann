# Hierarchical Boltzmann Simulations

This is a **very *basic proof-of-concept* research code** that solves the linear 
steady Boltzmann equation with Hermite-discretizations and an implicit
discontinuous Galerkin method in two dimensions on unstructured grids. 
The hierarchical nature allows to produce both accurate and very efficient 
solutions of classical fluid dynamics and precisely resolved Boltzmann simulations.

The code comes with **rather limited** documentation and user friendliness. For questions
please contact

    Manuel Torrilhon  - mt@mathcces.rwth-aachen.de
    Neeraj Sarna      - sarna@mathcces.rwth-aachen.de

## Requirements

The code has been tested on Unix- and OSX-systems using recent Gnu- and Intel-
compilers. The following additional libraries are required for compilation:
* C++-template library for linear algebra **Eigen**, see http://eigen.tuxfamily.org
* The Mesh-oriented Database **MOAB**, see http://sigma.mcs.anl.gov/moab-library

The linear system solution is based on the sparse-direct-solver pardiso which
can be used in two versions:
* Intel's math kernel library **Intel MKL**, see https://software.intel.com/en-us/node/470282
* original **Pardiso Solver Project**, see http://www.pardiso-project.org  

Some of the analytical solutions require Bessel functions provided by
* C++ libraries **Boost**, see http://www.boost.org

Some basic evaluation and visualization of the results is based on scripts
using
* **Mathematica**, see https://www.wolfram.com/mathematica 

## Basic Usage

1. Download the repository
2. Install the required libraries
3. Consult and modify the directory definitions in `MakefileMKL.def` or 
`MakefilePARDISO.def` depending on your installation.
4. Change into one of the application folders `System*`, the most general is 
`SystemGxy`, and compile with `make`
5. Check out the `Readme.md` of the application for instructions how to run 
the code

Alternatively, open the mathematica notebooks provided in the Evaluation folders
to explore precomputed results.

## Reference

This repository supplements the publication
*  M. Torrilhon and N. Sarna, 
   *Hierarchical Boltzmann Simulations and Model Error Estimation*,
   (2016), submitted

which also contains additional details of the methods used in the
code.

  