/*************************************************************************\
 SystemA.cpp  - implements setup of model system (A)
\*************************************************************************/
#include <stdio.h>
#include <iostream>
#include <math.h>
using namespace std;

#include "EigenSetup.h"
#include <Eigen/Eigenvalues>
#include "Mesh.h"
#include "System.h"
#include "SolutionA.h" 
#include "Tools.h"


// writes logfile message
string System::InfoText()
{
  return( strprintf( "System = %s, nEqn = %d, tau = %.2f, zeta = %.2f, kappa = %.2f, chi = %.2f, theta0 = %.2f, theta1 = %.2f, uW = %.2f, A0 = %.2f, A1 = %.2f, A2 = %.2f",
	                   SystemText.c_str(), nEqn, tau, zeta, kappa, chi, theta0, theta1, uW, A0, A1, A2 ) );  
};

// constructor with command line arguments
System::System( int argc, const char **argv )
{
  startTimer(">> System Setup \n"); 
  
//---------------------------------------------------------------------------------
// main system setup, maybe adjusted by user  
  tau = 1.0;
  if(exists_arg(argc,argv,"-tau")) tau = (double)get_float_arg(argc,argv,"-tau");
        
  kappa = 0.0;
  if(exists_arg(argc,argv,"-kappa")) kappa = (double)get_float_arg(argc,argv,"-kappa");
  
  outputFlag = 0;
  if(exists_arg(argc,argv,"-output")) outputFlag = (int)get_float_arg(argc,argv,"-output");  

  SystemText = "A(radial/BCodd)";

  setEQNData();

  zeta = 1.0; // legacy parameter in analytical solution
  chi = 1.0;
  
  if( SystemText.find("mixed") != -1 ) {
    theta0 = 1.0;
    theta1 = 0.5;
    uW = 5.0;  
    A0 = 2.0;
    A1 = 0.0;
    A2 = -1.0;
  };
  if( SystemText.find("radial0") != -1 ) {
    theta0 = 1.0;
    theta1 = 0.5;
    uW = 0.0;  
    A0 = 0.0;
    A1 = 0.0;
    A2 = 0.0;
  };
  if( SystemText.find("radial") != -1 ) {
    theta0 = 1.0;
    theta1 = 0.5;
    uW = 0.0;  
    A0 = 2.0;
    A1 = 0.0;
    A2 = -1.0;
  };
  BC.coeffRef(0,3) = kappa;  // kappa = 0 or 1
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

  cout << "System A \n";
  cout << "Number of equations " << nEqn << "\n";
  cout << "Number of boundary conditions " << nBC << "\n";
  cout << "maximal charac. velocity " << maxEV << "\n";

  setupExact();  
};

// gives a 'force' vector for the entire rhs of the system
MatrixXd System::Force( MatrixXd& pos )
{
  MatrixXd res(nEqn,1); 
  res = MatrixXd::Zero(nEqn,1);
  double r = pos.norm();
  
  res(0) = A0 + A2*r*r + A1*pos(0)/r;
  return( res );  
};

// sets the boundary data according to the boundary ID from the mesh file
void System::setBCData( MatrixXd& pos, MatrixXd& normal, int boundaryID )
{
  MatrixXd g(nBC,1);
  double norm = pos.norm(), chiW, vtW, thetaW;
  
  switch(boundaryID) {
  case 10:// inner circle
    chiW = chi;
    thetaW = theta0;
    vtW = uW*normal(1);
    break;
  case 20:// outer circle
    chiW = chi;
    thetaW = theta1;
    vtW = 0.0;
    break;
  }; 
 
  XBC = chiW*BC;
  g = -vtW*XBC.col(2) -thetaW*XBC.col(0);  
  for( int ix=0; ix<OddVar.size(); ix++ ) XBC.coeffRef(ix,OddVar(ix)) = BC.coeffRef(ix,OddVar(ix));

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
  string Text, text1, text2;
   
  text1 = ComputeError( Mesh );
  Text = text1;

  if( outputFlag == 1 ) {
    text2 = FileOutput( Mesh );
    Text = strprintf( "%s\n%s", text1.c_str(), text2.c_str() );
  };

  return( Text );
};
  
// compute L2/Linf error of temperature and heat-flux
string System::ComputeError( MeshCore *Mesh )
{  
  int length;
  Tag Solution;
  double err_theta_inf = 0.0, err_theta_L2 = 0.0, err_sx_inf = 0.0, err_sx_L2 = 0.0;
  double theta_max = 0.0, sx_max = 0.0;

  startTimer(" compute error \n"); 
  Mesh->tag_get_handle("Solution Elements",Solution);
  Mesh->tag_get_length(Solution,length);
  
  for (Range::iterator elem = Mesh->Elements.begin(); elem != Mesh->Elements.end(); elem++) { 
    const EntityHandle *vertices;
    int type;	
    Mesh->get_connectivity(*elem, vertices, type, true);
    MatrixXd vertP(3,type);
    Mesh->get_coords(vertices,type,&vertP(0));
    MatrixXd pos = vertP.rowwise().mean();    
    MatrixXd solution(length,1);
    Mesh->tag_get_data(Solution,&(*elem),1,&solution(0)); 
    
    double theta = solution(0), sx = solution(1);
    double err_theta = fabs(theta - solX->theta(pos(0),pos(1)));
    double err_sx = fabs(sx - solX->s_x(pos(0),pos(1)));
    
    err_theta_inf = max(err_theta_inf,err_theta);
    err_sx_inf = max(err_sx_inf,err_sx);
    err_theta_L2 += err_theta*err_theta;
    err_sx_L2 += err_sx*err_sx;
    
    theta_max = max(theta_max,fabs(solX->theta(pos(0),pos(1))));
    sx_max = max(sx_max,fabs(solX->s_x(pos(0),pos(1))));
  };
  
  MatrixXd err(4,1);
  err << err_theta_inf/theta_max, Mesh->hMax*sqrt(err_theta_L2)/theta_max, err_sx_inf/sx_max, Mesh->hMax*sqrt(err_sx_L2)/sx_max;
  stopTimer();  
  return( strprintf( "Errors = theta/sx, theta-Linf = %f, theta-L2 = %f, sx-Linf = %f, sx-L2 = %f", err(0), err(1), err(2), err(3) ) );
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
    
    double theta = solution(0), qx = solution(1), qy = solution(2);
     
    fprintf( fp, "%20.8f", pos(0) );
    fprintf( fp, "%20.8f", pos(1) );
    fprintf( fp, "%20.8f", theta );
    fprintf( fp, "%20.8f", qx );
    fprintf( fp, "%20.8f", qy );
    fprintf( fp, "\n" );          
  };
  
  fclose( fp );
  stopTimer();
  return( strprintf( "ResultFile = '%s'", filename.c_str() ) );
};


