/*************************************************************************\
 pardisoSolvePARDISO.cpp  - interface routine to original direct linear
                            solver of Pardiso.
This is based on the package available at 
http://www.pardiso-project.org/          
\*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"


// pardiso prototypes
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                            double *, int    *,    int *, int *,   int *, int *,
                            int *, double *, double *, int *, double *);

// conducts a solution phase of the direct solver pardiso
void PardisoSolve(INT mtype, INT n, INT *ia, INT *ja, double *A, int nz,
									double *b, double *x, void *pt, INT phase, INT *iparm, double *dparm )
{      
  int maxfct = 1;  // Maximum number of numerical factorizations. 
  int mnum   = 1;  // Which factorization to use.    
  int msglvl = 0;  // = 1 Print statistical information 
  int error  = 0;  // Initialize error flag 
  int solver = 0;  // use sparse direct solver 
  int nrhs = 1;    // Number of right hand sides
  double   ddum;   // Double dummy 
  int      idum;   // Integer dummy

  printf( " PardisoSolve :" );
  switch( phase ) 
  {
    case( -11 ):    // Pardiso control parameters.
      printf( " Setup Phase" );

      pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);
      printf("\n");

      if (error != 0) {
          if (error == -10 ) printf("No license file found \n");
          if (error == -11 ) printf("License is expired \n");
          if (error == -12 ) printf("Wrong username or hostname \n");
          exit(0);
      }   

      int num_procs; // Numbers of processors, value of OMP_NUM_THREADS 
      char *var;
      var = getenv("OMP_NUM_THREADS");
      if(var != NULL)
          sscanf( var, "%d", &num_procs );
      else {
          printf("Set environment OMP_NUM_THREADS to 1");
          exit(1);
      }
      iparm[2]  = num_procs;
      iparm[7] = 200; // Max numbers of iterative refinement steps  
      iparm[9] = 10; // Perturb the pivot elements with 1E-8
      iparm[12] = 0; // 1 = Maximum weighted matching algorithm is switched-on (default for non-symmetric)  

      //Convert matrix from 0-based C-notation to Fortran 1-based notation.         
      for(int i = 0; i < n+1; i++) ia[i] += 1; 
      for(int i = 0; i < nz; i++) ja[i] += 1;

    break;
    case( 11 ):    // Symbolic Factorization
      printf( " Symbolic Factorization \n" );
	    pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
        &n, A, (INT *)ia, (INT *)ja, &idum, &nrhs,
		    iparm, &msglvl, &ddum, &ddum, &error, dparm);
	    if (error != 0) {
		    printf("\nERROR during symbolic factorization: %d\n", error);
		    exit(1);
	    }
    break;
    case( 22 ):    // Numerical Factorization
      printf( " Numerical Factorization \n" );
	    pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
        &n, A, (INT *)ia, (INT *)ja, &idum, &nrhs,
		    iparm, &msglvl, &ddum, &ddum, &error, dparm);
	    if (error != 0) {
		    printf("\nERROR during numerical factorization: %d\n", error);
		    exit(2);
	    }
    break;
    case( 33 ):   //Back substitution and iterative refinement
      printf( " Back Substitution and Iterative Refinement \n" );
	    pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
        &n, A, (INT *)ia, (INT *)ja, &idum, &nrhs,
		    iparm, &msglvl, b, x, &error, dparm);
	    if (error != 0) {
		    printf("\nERROR during solution: %d\n", error);
		    exit(3);
	    }
    break;
    case( -1 ): // Release internal memory. 
      printf( " Release Internal Memory \n" );
	    pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
        &n, A, (INT *)ia, (INT *)ja, &idum, &nrhs,
		    iparm, &msglvl, &ddum, &ddum, &error, dparm);

    //Convert matrix back to 0-based C-notation. 
    for(int i = 0; i < n+1; i++) ia[i] -= 1;
    for(int i = 0; i < nz; i++) ja[i] -= 1;
		break;
	};
};


