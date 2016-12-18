/*************************************************************************\
 Equations.cpp  - equation system matrices for linear Navier-Stokes-
                  Fourier system. 
\*************************************************************************/
#include <stdio.h>
#include <iostream>
#include <math.h>
using namespace std;

#include "EigenSetup.h"
#include <Eigen/Eigenvalues>
#include "Mesh.h"
#include "System.h"
#include "Tools.h"

// sets the matrices for the equations and boundary conditions 
void System::setEQNData()
{
  VectorXi Ai, Aj, Pi, Pj, BCi, BCj;
  MatrixXd Aval, Pval, BCval;
  
  // p, (vx, vy), theta, (sig_xx, sig_xy, sig_yy), (qx, qy)
  nEqn = 9;
  nBC = 3;

  Ai.resize(11);
  Aj.resize(11);
  Aval.resize(11,1);
  Ai <<     0,   1,   1,   1,   2,   3,   3,     4,   5,      6,  7;
  Aj <<     1,   4,   3,   0,   5,   7,   1,     1,   2,      1,  3;
  Aval << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 4.0/3, 1.0, -2.0/3, 2.5;
   
  Pi.resize(5);
  Pj.resize(5);
  Pval.resize(5,1);
  Pi << 4,5,6,7,8;
  Pj << 4,5,6,7,8;
  Pval << 1,1,1,2.0/3,2.0/3;
  
  BCi.resize(7);
  BCj.resize(7);
  BCval.resize(7,1);
  BCi << 0,0,0,1,1,2,2;
  BCj << 0,1,4,2,5,3,7;
  BCval << 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0;
  OddVar.resize(nBC);
  OddVar << 1,5,7; 
  TensorDegree.resize(5);
  TensorDegree << 0,1,0,2,1; 
  
  SpInit( nEqn, nEqn, Ai, Aj, Aval, Ax );
  SpInit( nEqn, nEqn, Pi, Pj, Pval, P );
  SpInit( nBC, nEqn, BCi, BCj, BCval, BC );
};

// generates the projector of normal coordinates 
SpMatrix System::Projector( const MatrixXd& normal )
{
  double nx = normal(0), ny = normal(1);
  int idx;
  MatrixXd P[8];
  P[0].resize(1,1);
  P[0] << 1.0;
  P[1].resize(2,2);
  P[1] << nx, ny, -ny, nx;
  P[2].resize(3,3);
  double n20 = nx*nx, n02 = ny*ny, n11 = nx*ny;
  P[2] << n20,2*n11,n02,-n11,-n02 + n20,n11,n02,-2*n11,n20;
  
  SpMatrix T(nEqn,nEqn);
  VectorXi ColSize(nEqn);
  int ix = 0;
  for( int i=0; i<TensorDegree.size(); i++ ) 
    for( int j=0; j<TensorDegree(i)+1; j++ )
      ColSize(ix++) = TensorDegree(i)+1;
      
  T.reserve(ColSize);
  ix = 0;
  for( int i=0; i<TensorDegree.size(); i++ ) {
    SpBlock( ix, P[TensorDegree(i)], T );
    ix += TensorDegree(i)+1;
  };
  
  T.makeCompressed();
  return( T );
};

// defines the inverse projector
SpMatrix System::invProjector( const MatrixXd& normal )
{
  MatrixXd mirrow(3,1);
  mirrow << normal(0,0),-normal(1,0),0;
  return( Projector( mirrow ) );
};

