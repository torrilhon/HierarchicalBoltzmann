/*************************************************************************\
 Equations.cpp  - equation system matrices for system (B) 
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
  
  // theta, (qx, qy), (R_xx, R_xy, R_yy), (psi_xxx,psi_xxy,psi_xyy,psi_yyy)
  nEqn = 10;
  nBC = 4;

  Ai.resize(15);
  Aj.resize(15);
  Aval.resize(15,1);
  Ai << 0,1,1,2,3,3,4,4,5,5,6,7,8,8,9;
  Aj << 1,3,0,4,6,1,7,2,1,8,3,4,5,3,4;
  Aval << 1.0,1.0,1.0,1.0,1.0,2.0/3,1.0,0.5,-1.0/3,1.0,6.0/5,16.0/15,2.0/3,-4.0/15,-4.0/5;
   
  Pi.resize(9);
  Pj.resize(9);
  Pval.resize(9,1);
  Pi << 1,2,3,4,5,6,7,8,9;
  Pj << 1,2,3,4,5,6,7,8,9;
  Pval << 1,1,1,1,1,1,1,1,1;
  
  BCi.resize(10);
  BCj.resize(10);
  BCval.resize(10,1);
  BCi << 0,0,0,1,1,1,2,2,3,3;
  BCj << 0,1,3,2,4,7,3,6,5,8;
  BCval << 1.0,-1.0,1.0,1.0,-1.0,2.0,1.0,-1.0,1.0,-1.0;
  OddVar.resize(nBC);
  OddVar << 1,4,6,8; 
  TensorDegree.resize(4);
  TensorDegree << 0,1,2,3; 
  
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

