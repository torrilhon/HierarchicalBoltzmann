/*************************************************************************\
 Equations.cpp  - equation system matrices for linear R13 system. 
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
  
  // rho, vx, vy, theta, sigma_xx, sigma_xy, sigma_yy, qx, qy, mxxx, mxxy, mxyy, myyy, Rxx, Rxy, Ryy
  nEqn = 16;
  nBC = 6;

  Ai.resize(29);
  Aj.resize(29);
  Aval.resize(29,1);
  Ai << 0, 1, 1, 1, 2, 3, 3, 4, 4, 4, 5, 5,  5,  6, 6, 6, 7, 7,  7, 8, 8, 9, 10, 11, 11, 12, 13, 14, 15;
  Aj << 1, 3, 0, 4, 5, 1, 7, 7, 1, 9, 8, 10, 2, 11, 1, 7, 4, 3, 13, 5, 14, 4, 5, 4,  6,  5,  7,  8,  7;
  Aval << 1., 1., 1., 1., 1., 1., 1., 8.0/15, 4.0/3, 1., 0.4, 1., 1., 1., -2.0/3, -4.0/15, 1., 2.5, 0.5, 
          1., 0.5, 1.8, 1.6, -0.4, 1.0, -1.2, 56.0/15, 2.8, -28.0/15;
   
  Pi.resize(12);
  Pj.resize(12);
  Pval.resize(12,1);
  Pi << 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
  Pj << 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
  Pval << 1., 1., 1., 2.0/3, 2.0/3, 1.5, 1.5, 1.5, 1.5, 7.0/6, 7.0/6, 7.0/6;
  
  BCi.resize(24);
  BCj.resize(24);
  BCval.resize(24,1);
  BCi << 0,0,0,  1,1,1, 1,  2,2,2, 2,  3,3,3, 3,  4,4,4, 4, 4,  5,5, 5, 5;
  BCj << 0,1,3,  2,5,8,10,  3,4,7,13,  3,4,9,13,  3,4,6,11,13,  2,8,10,14;
  BCval << 1.0,-1.0,1.0,
           1.0,-1.0,0.2,1.0,
           2.0,0.5,-1.0,0.4,
           -0.4,1.4,-1.0,-0.08,
           0.2,-0.2,1.0,-1.0,0.04,
           -1.0,2.2,-1.0,-1.0;
  OddVar.resize(nBC);
  OddVar << 1,5,7,9,11,14; 
  TensorDegree.resize(7);
  TensorDegree << 0,1,0,2,1,3,2; 
  
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
  P[3].resize(4,4);
  double n30 = nx*n20, n21 = ny*n20, n12 = ny*n11, n03 = ny*n02;
  P[3] << n30,3*n21,3*n12,n03,-n21,-2*n12+n30,-n03+2*n21,n12,n12,n03-2*n21,-2*n12+n30,n21,-n03,3*n12,-3*n21,n30;
  
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

