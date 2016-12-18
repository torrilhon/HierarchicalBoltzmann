/*************************************************************************\
 Equations.cpp  - equation system matrices for system (A) 
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
  
  // theta, (qx, qy), (R_xx, R_xy, R_yy)
  nEqn = 6;
  nBC = 2;

  Ai.resize(7);
  Aj.resize(7);
  Aval.resize(7,1);
  Ai << 0,1,1,2,3,4,5;
  Aj << 1,0,3,4,1,2,1;
  Aval << 1.0,1.0,1.0,1.0,2.0/3,0.5,-1.0/3;
  
  Pi.resize(5);
  Pj.resize(5);
  Pval.resize(5,1);
  Pi << 1,2,3,4,5;
  Pj << 1,2,3,4,5;
  Pval << 1,1,1,1,1;
  
  BCi.resize(5);
  BCj.resize(5);
  BCval.resize(5,1);
  BCi << 0,0,0,1,1;
  BCj << 0,1,3,2,4;
  BCval << 1.0,-1.0,1.0,1.0,-1.0;

  OddVar.resize(nBC);
  OddVar << 1,4; 
  TensorDegree.resize(3);
  TensorDegree << 0,1,2; 
  
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

