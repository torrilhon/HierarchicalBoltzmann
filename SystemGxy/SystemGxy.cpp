/*************************************************************************\
 SystemGxy.cpp  - implements setup and output for hierarchical Boltzmann  
                  discretizations.
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


// writes logfile message
string System::InfoText()
{
  return(strprintf("System = %s, nMoments = %d, nEqn = %d, tau = %.3f, chi = %.3f, theta0 = %.3f, theta1 = %.3f, v0 = %.3f, p0 = %.3f, eps = %.2e",
	                 SystemText.c_str(), nMoments, nEqn, tau, chi, theta0, theta1, v0, p0, eps ));
};

// constructor with command line arguments
System::System( int argc, const char **argv )
{
  startTimer(">> System Setup \n"); 

//---------------------------------------------------------------------------------
// main system setup, maybe adjusted by user  
  tau = 1.0;
  if(exists_arg(argc,argv,"-tau")) tau = (double)get_float_arg(argc,argv,"-tau");

  nMoments = 20;
  if(exists_arg(argc,argv,"-moments")) nMoments = (int)get_float_arg(argc,argv,"-moments");
  
  outputFlag = 1;
  if(exists_arg(argc,argv,"-output")) outputFlag = (int)get_float_arg(argc,argv,"-output");  

  SystemText = "Gxy(heat+flow/BCodd)";

  setEQNData( nMoments );

  chi = 1.0;    // accommodation coefficient 
  theta0 = 1.0; // inner wall temperature 
  theta1 = 1.0; // outer wall tempearture
  p0 = 1.0;     // inflow pressure
  v0 = 1.0;     // inflow velocity (not used)
  eps = 1e-5;   // boundary regularization for velocity/pressure
//---------------------------------------------------------------------------------
  
// system preparation for numerical method  
  P = 1.0/tau*P;

  Ay.resize(nEqn,nEqn);
  Matrix<double,2,1> normal(0.0,1.0);
  Ay = invProjector(normal) * Ax * Projector(normal);
  Ay.prune(1,1e-13);
  stopTimer();
  
  startTimer(" characteristic setup \n"); 
  EigenSolver<MatrixXd> ES(Ax);
  MatrixXd vecs = ES.pseudoEigenvectors();
  VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal();
  double maxEV = vals.cwiseAbs().maxCoeff();
  Aminus1D = vecs*(vals.cwiseAbs()-vals).asDiagonal()*vecs.inverse();

  X0.resize(nEqn,nBC);
  for( int i=0; i<nBC; i++ ) X0.col(i) = VectorXd::Unit(nEqn,OddVar(i)); // odd BCs
  
  Ordering Order(vals);
  //for( int i=0; i<nBC; i++ ) X0.col(i) = vecs.col(Order.index(i)); // Riemann BCs
  
  XBC.resize(nEqn,nEqn);
  BCrhs.resize(nEqn,1); 
  stopTimer();

  cout << "Moment theory " << nMoments << "\n";
  cout << "Number of equations " << nEqn << "\n";
  cout << "Number of boundary conditions " << nBC << "\n";
  cout << "maximal charac. velocity " << maxEV << "\n";
  
};

// gives a 'force' vector for the entire rhs of the system
MatrixXd System::Force( MatrixXd& pos )
{
  MatrixXd rhs(nEqn,1); 
  rhs = MatrixXd::Zero(nEqn,1);
  double r = pos.norm();
  
  return( rhs );  
};

// sets the boundary data according to the boundary ID from the mesh file
void System::setBCData( MatrixXd& pos, MatrixXd& normal, int boundaryID )
{
  MatrixXd g(nBC,1);
  double norm = pos.norm(), chiW, eps_p, vtW, vnW, pW, thetaW;
  
  switch(boundaryID) {
  case 10:// inflow
    chiW = 100*chi;
    eps_p = eps;
    vnW = -v0;
    vtW = 0.0;
    pW = p0;
    thetaW = theta0;
    break;
  case 20:// outflow
    chiW = 100*chi;
    eps_p = eps;
    vnW = v0;
    vtW = 0.0;
    pW = -p0;
    thetaW = theta0;
    break;
  case 30:// wall
    chiW = chi;
    eps_p = 1.0/eps;
    vnW = 0.0;
    vtW = 0.0;
    pW = 0.0;
    thetaW = theta0;
    break;
  case 40:// inner walls
  case 50:
  case 60:
    chiW = chi;
    eps_p = 1.0/eps;
    vnW = 0.0;
    vtW = 0.0;
    pW = 0.0;
    thetaW = theta1;
    break;
  }; 
 
  XBC = chiW*BC;
  g = -vtW*XBC.col(2) -(-3.0/2)*thetaW*XBC.col(3);
  for( int ix=0; ix<OddVar.size(); ix++ ) XBC.coeffRef(ix,OddVar(ix)) = BC.coeffRef(ix,OddVar(ix));
  XBC.coeffRef(0,1) = -eps_p*BC.coeffRef(0,1);
  g(0) = +vnW*eps_p*BC.coeffRef(0,1) -pW*chiW*BC.coeffRef(0,0);

  BCrhs = invProjector(normal)*X0*(XBC*X0).inverse()*g;
  XBC = invProjector(normal)*X0*(XBC*X0).inverse()*XBC*Projector(normal);
};
 

// upwind flux matrix
MatrixXd System::Aminus( const MatrixXd& normal )
{
  return(invProjector(normal)*Aminus1D*Projector(normal)  );
};

// post processing the solution
string System::PostProcess( MeshCore *Mesh )
{
  if( outputFlag == 1 ) 
    return( FileOutput( Mesh ) );

  return( string("") );
};

// simple data output into a time-stamped dat-file
string System::FileOutput( MeshCore *Mesh )
{
  string filename;
  int length;
  FILE *fp;
  Tag Solution;

  filename = strprintf( "../output/result%s.dat", TimeStamp().c_str() );

  startTimer( filename );
  cout << endl;
  fp = fopen( filename.c_str(), "w" );

  Mesh->tag_get_handle("Solution Vertices",Solution);
  Mesh->tag_get_length(Solution,length);
  
  for (Range::iterator vert = Mesh->Vertices.begin(); vert != Mesh->Vertices.end(); vert++) {
    MatrixXd pos(3,1), solution(length,1);
    Mesh->get_coords(&(*vert),1,&pos(0));
    Mesh->tag_get_data(Solution,&(*vert),1,&solution(0)); 
    
    double rho = solution(0), vx = solution(1), vy = solution(2), theta = (-2.0/3)*solution(3), 
           sigxx = solution(4), sigxy = solution(5), sigyy = solution(6), 
           qx = (-1)*solution(7), qy = (-1)*solution(8);
    double stress = sqrt(3*(sigxx*sigxx + sigxy*sigxy + sigxx*sigyy + sigyy*sigyy));  
     
    fprintf( fp, "%20.8f", pos(0) );
    fprintf( fp, "%20.8f", pos(1) );
    fprintf( fp, "%20.8f", rho );
    fprintf( fp, "%20.8f", vx );
    fprintf( fp, "%20.8f", vy );
    fprintf( fp, "%20.8f", theta );
    fprintf( fp, "%20.8f", stress );
    fprintf( fp, "%20.8f", qx );
    fprintf( fp, "%20.8f", qy );
    fprintf( fp, "\n" );          
  };
  
  fclose( fp );
  stopTimer();
  return( strprintf( "ResultFile = '%s'", filename.c_str() ) );
};


