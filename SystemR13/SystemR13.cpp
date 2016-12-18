/*************************************************************************\
 SystemR13.cpp  - implements setup of R13 system.
\*************************************************************************/
#include <stdio.h>
#include <iostream>
#include <math.h>
using namespace std;

#include "EigenSetup.h"
#include <Eigen/Eigenvalues>
#include "Mesh.h"
#include "System.h"
#include "SolutionR13.h" 
#include "Tools.h"


// writes logfile message
string System::InfoText()
{
  return(strprintf( "System = %s, nEqn = %d, tau = %.2f, chi = %.2f, theta0 = %.2f, theta1 = %.2f, v0 = %.2f, A0 = %.2f, A1 = %.2f, A2 = %.2f, eps = %.2e",
	                  SystemText.c_str(), nEqn, tau, chi, theta0, theta1, v0, A0, A1, A2, eps ) );
};

// constructor with command line arguments
System::System( int argc, const char **argv )
{
  startTimer(">> System Setup \n"); 
  
//---------------------------------------------------------------------------------
// main system setup, maybe adjusted by user  
  tau = 1.0;
  if(exists_arg(argc,argv,"-tau")) tau = (double)get_float_arg(argc,argv,"-tau");
        
  outputFlag = 0;
  if(exists_arg(argc,argv,"-output")) outputFlag = (int)get_float_arg(argc,argv,"-output");  

  SystemText = "R13(heat+flow/BCodd)";
  
  setEQNData();

  chi = 1.0;    // accommodation coefficient 
  theta0 = 1.0; // inner wall temperature 
  theta1 = 2.0; // outer wall tempearture
  v0 = 1.0;     // inflow velocity (not used)
  p0 = 0.27;    // pressure constant
  A0 = 0.0;     // not supported (legacy)
  A1 = 0.0;     // not supported (legacy)
  A2 = 0.0;     // not supported (legacy)
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

  cout << "R13 System\n";
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
  
  res(3) = A0 + A2*r*r + A1*pos(0)/r;
  return( res );  
};

// sets the boundary data according to the boundary ID from the mesh file
void System::setBCData( MatrixXd& pos, MatrixXd& normal, int boundaryID )
{
  MatrixXd g(nBC,1);
  double norm = pos.norm(), chiW, eps_p, eps_n, vtW, vnW, pW, thetaW;
  
  switch(boundaryID) {
  case 10:// inner circle
    chiW = chi;
    //eps_n = eps;
    eps_p = 1.0/eps;
    vnW = 0.0;
    vtW = 0.0;
    pW = 0.0;
    thetaW = theta0;
    break;
  case 20:// outer circle
    chiW = chi;
    //eps_n = 1.0;
    eps_p = eps;
    vnW = v0*normal(0);
    vtW = -v0*normal(1);
    pW = -p0*normal(0);
    thetaW = theta1;
    break;
  }; 
 
  XBC = chiW*BC;
  g << vnW*eps_p - pW*chiW, -vtW*chiW, -2*thetaW*chiW, 0.4*thetaW*chiW, -0.2*thetaW*chiW, vtW*chiW;
  for( int ix=0; ix<OddVar.size(); ix++ ) XBC.coeffRef(ix,OddVar(ix)) = BC.coeffRef(ix,OddVar(ix));
  XBC.coeffRef(0,1) = -eps_p;
  
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
  double err_ux_inf = 0.0, err_ux_L2 = 0.0, err_sigxy_inf = 0.0, err_sigxy_L2 = 0.0;
  double theta_max = 0.0, sx_max = 0.0, ux_max = 0.0, sigxy_max = 0.0;

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
    
    double theta = solution(3), ux = solution(1), sigxy = solution(5), sx = solution(7);
    double err_theta = fabs(theta - solX->theta(pos(0),pos(1)));
    double err_ux = fabs(ux - solX->u_x(pos(0),pos(1)));
    double err_sx = fabs(sx - solX->s_x(pos(0),pos(1)));
    double err_sigxy = fabs(sigxy - solX->sig_xy(pos(0),pos(1)));
    
    err_theta_inf = max(err_theta_inf,err_theta);
    err_sx_inf = max(err_sx_inf,err_sx);
    err_ux_inf = max(err_ux_inf,err_ux);
    err_sigxy_inf = max(err_sigxy_inf,err_sigxy);
    err_theta_L2 += err_theta*err_theta;
    err_sx_L2 += err_sx*err_sx;
    err_ux_L2 += err_ux*err_ux;
    err_sigxy_L2 += err_sigxy*err_sigxy;
    
    theta_max = max(theta_max,fabs(theta));
    sx_max = max(sx_max,fabs(sx));
    ux_max = max(ux_max,fabs(ux));
    sigxy_max = max(sigxy_max,fabs(sigxy));
  };
  
  MatrixXd err(8,1);
  err << err_theta_inf/theta_max, Mesh->hMax*sqrt(err_theta_L2)/theta_max, err_sx_inf/sx_max, Mesh->hMax*sqrt(err_sx_L2)/sx_max,
         err_ux_inf/ux_max, Mesh->hMax*sqrt(err_ux_L2)/ux_max, err_sigxy_inf/sigxy_max, Mesh->hMax*sqrt(err_sigxy_L2)/sigxy_max;
  stopTimer();  
  return(strprintf("Errors = theta/sx/ux/sigxy, theta-Linf = %f, theta-L2 = %f, sx-Linf = %f, sx-L2 = %f, ux-Linf = %f, ux-L2 = %f, sigxy-Linf = %f, sigxy-L2 = %f", 
                   err(0), err(1), err(2), err(3), err(4), err(5), err(6), err(7) ) );
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
    
    double rho = solution(0), vx = solution(1), vy = solution(2), theta = solution(3), 
           sigxx = solution(4), sigxy = solution(5), sigyy = solution(6), 
           qx = solution(7), qy = solution(8);
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

