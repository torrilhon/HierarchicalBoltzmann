/*************************************************************************\
 pardisoSolveMKL.cpp  - interface routine to the INtel-MKL version of 
                        the direct linear solver of Pardiso.
This is based on the Intel Math-Kernel-Library extension 
https://software.intel.com/en-us/node/470282          
\*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mkl.h>
#include "defs.h"

// conducts a solution phase of the direct solver pardiso
void PardisoSolve(INT mtype, INT n, INT *ia, INT *ja, double *A, int nz,
                  double *b, double *x, void *pt, INT phase, INT *iparm, double *dparm)
{  
	INT i, nrhs;
	INT maxfct, mnum, msglvl, error;       
	INT idum;    // Integer dummy 
	double ddum;  // Double dummy 
	nrhs = 1;    // number of rhs
	maxfct = 1;  // Maximum number of numerical factorizations
	mnum = 1;    // Which factorization to use.
	msglvl = 0;  // 1 = Print statistical information in file 
	error = 0;   // Initialize error flag 
    
  printf( " PardisoSolve :" );
  switch( phase ) 
  {
    case( -11 ):    // Pardiso control parameters.
      printf( " Setup Phase \n" );     
	      
    	for (i = 0; i < 64; i++) iparm[i] = 0;             

	    iparm[0] = 1; // No solver default  
	    iparm[1] = 3; // Fill-in reordering from METIS  
	    iparm[2] = 1; // Numbers of processors, value of OMP_NUM_THREADS  
	    iparm[3] = 0; // No iterative-direct algorithm  
	    iparm[4] = 0; // No user fill-in reducing permutation  
	    iparm[5] = 0; // 0 = Write solution into x  
	    iparm[6] = 0; // Not in use  
	    iparm[7] = 200; // Max numbers of iterative refinement steps  
	    iparm[8] = 0; // Not in use  
	    iparm[9] = 8; // Perturb the pivot elements with 1E-8  
	    iparm[10] = 1; // 1 = Use nonsymmetric permutation and scaling MPS  
	    iparm[11] = 0; // Not in use  
	    iparm[12] = 0; // 1 = Maximum weighted matching algorithm is switched-on (default for non-symmetric)  
	    iparm[13] = 0; // Output: Number of perturbed pivots  
	    iparm[14] = 0; // Not in use  
	    iparm[15] = 0; // Not in use  
	    iparm[16] = 0; // Not in use  
	    iparm[17] = -1; // Output: Number of nonzeros in the factor LU  
	    iparm[18] = -1; // Output: Mflops for LU factorization  
	    iparm[19] = 0; // Output: Numbers of CG Iterations  
	    iparm[34] = 1; // 0= one-based indices, 1= zero-besed indices 
 	    iparm[59] = 0; // 0= in-core, 2= out-of-core  

    break;
    case( 11 ):    // Symbolic Factorization
      printf( " Symbolic Factorization \n" );
	    pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
        &n, A, (INT *)ia, (INT *)ja, &idum, &nrhs,
		    iparm, &msglvl, &ddum, &ddum, &error);
	    if (error != 0) {
		    printf("\nERROR during symbolic factorization: %d\n", error);
		    exit(1);
	    }
    break;
    case( 22 ):    // Numerical Factorization
      printf( " Numerical Factorization \n" );
	    pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
        &n, A, (INT *)ia, (INT *)ja, &idum, &nrhs,
		    iparm, &msglvl, &ddum, &ddum, &error);
	    if (error != 0) {
		    printf("\nERROR during numerical factorization: %d\n", error);
		    exit(2);
	    }
    break;
    case( 33 ):   //Back substitution and iterative refinement
      printf( " Back Substitution and Iterative Refinement \n" );
	    pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
        &n, A, (INT *)ia, (INT *)ja, &idum, &nrhs,
		    iparm, &msglvl, b, x, &error);
	    if (error != 0) {
		    printf("\nERROR during solution: %d\n", error);
		    exit(3);
	    }
    break;
    case( -1 ): // Release internal memory. 
      printf( " Release Internal Memory \n" );
	    pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
        &n, A, (INT *)ia, (INT *)ja, &idum, &nrhs,
		    iparm, &msglvl, &ddum, &ddum, &error);
		break;
	};
};
