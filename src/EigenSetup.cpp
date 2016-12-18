/*************************************************************************\
 EigenSetup.cpp  - defines some helpfull initialization routines. 
 This is based on the Eigen template library, see
 http://eigen.tuxfamily.org/
\*************************************************************************/
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
using namespace std;

#include "EigenSetup.h"

// sparse matrix initialization from triplet arrays
void SpInit( int n, int m, const VectorXi& Ai, const VectorXi& Aj, const MatrixXd& Aval, SpMatrix& A )
{
  VectorXi ColSize(n);
  A.resize(n,m);
  ColSize.setZero();
  for( int ix=0; ix < Ai.size(); ix++ ) ColSize(Ai(ix))++; 
  A.reserve(ColSize);
  for( int ix=0; ix < Ai.size(); ix++ ) A.coeffRef(Ai(ix),Aj(ix)) = Aval(ix,0);
  A.prune(1,1e-13);  
};

// block matrix initialization
void SpBlock( int idx, const MatrixXd& P, SpMatrix& T )
{
  for( int i=0; i < P.rows(); i++ )
    for( int j=0; j < P.cols(); j++ )
      T.coeffRef(idx+i,idx+j) = P(i,j); 
};

