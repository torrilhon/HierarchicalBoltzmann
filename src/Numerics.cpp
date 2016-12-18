/*************************************************************************\
 Numerics.cpp  - assembles the DG matrix based on the provided system 
                 and solves the linear equations.
\*************************************************************************/
#include <stdio.h>
#include <iostream>
#include <math.h>
using namespace std;

#include "EigenSetup.h"
#include "Mesh.h"
#include "System.h"
#include "Numerics.h"
#include "Tools.h"

// declaration of the generic pardiso solver routine (see pardisoSolveXYZ.cpp)
void PardisoSolve(INT mtype, INT n, INT *ia, INT *ja, double *A, int nz,
                  double *b, double *x, void *pt, INT phase, INT *iparm, double *dparm);


// writes logfile message
string Numerics::InfoText()
{
  return( strprintf( "Method = DG, order = %d, nBasis = %d, nTot = %d, nnz = %ld, nnzRatio = %.2f, res = %.2e, nIter = %d", 
                     order, nBasis, nTot, nnz, 1.0*nnz/nTot, res, nIter ) );
};

// initialization and memory allocation
Numerics::Numerics( MeshCore *mesh, System *sys, int argc, const char **argv )
{
  startTimer(">> Numerics Setup\n");
  Sys = sys;
  Mesh = mesh;  
 
  order = 2;
  if(exists_arg(argc,argv,"-order")) order = (int)get_float_arg(argc,argv,"-order");
 
  if( Mesh->Type == 3 ) nBasis = (order+1)*(order+2)/2;
  if( Mesh->Type == 4 ) nBasis = (order+1)*(order+1);
  nTot = Mesh->Elements.size()*Sys->nEqn*nBasis;
  
  nnz = 0;
  res = 0.0;
  nIter = -1;
  stopTimer();
  
  startTimer(" allocate stencil size\n");
  VectorXi StencilSize(nTot);// the knowledge of the sparsity pattern accelerates the matrix build-up
  StencilSize.setZero(nTot);
  StencilSize = VectorXi::Constant(nTot,4*Sys->nEqn*nBasis);// precise band width needs to be improved!!
  stopTimer();
  
  startTimer(" allocate sparse matrices \n");
  Atotal.resize(nTot,nTot);
  Atotal.reserve(StencilSize);
  sol.resize(nTot,1);
  sol.setZero(nTot,1);
  rhs.resize(nTot,1);
  rhs.setZero(nTot,1);
  stopTimer();  
  
  cout << "Order of Elements " << order << "\n";
  cout << "Number of basis functions " << nBasis << "\n";
  cout << "Total number of unknowns " << nTot << "\n"; 
}

// core routine of DG method which assembles the system matrix into sparse format
void Numerics::Assemblation()
{
  startTimer(">> Assemblation - element-wise\n"); 
  
  // global element loop
  for (Range::iterator elem = Mesh->Elements.begin(); elem != Mesh->Elements.end(); elem++) { 
    std::vector<EntityHandle> vertices;
    int type;
    Mesh->get_connectivity(&(*elem), 1, vertices, true );
    type = vertices.size();
    vertices.resize(type+1);
    vertices[type] = vertices[0];
    
    MatrixXd pos(3,nBasis);
    Mesh->get_coords(&(vertices[0]),type,&pos(0));
    LagrangePts( pos, order, type );
    
    MatrixXd AssemMatX = AssembleDiffX(pos,order,type);
    MatrixXd AssemMatY = AssembleDiffY(pos,order,type);
    MatrixXd AssemMass = AssembleMass(pos,order,type);
    
    int offset = Mesh->get_ID(*elem)*Sys->nEqn*nBasis;
    
    // volume integrals
    for( int i=0; i<nBasis; i++ ) {
      for( int j=0; j<nBasis; j++ ) {

        for( int k=0; k<Sys->Ax.outerSize(); ++k)
          for (SpMatrix::InnerIterator it(Sys->Ax,k); it; ++it ) 
            Atotal.coeffRef(offset+it.row()*nBasis+i,offset+it.col()*nBasis+j) += it.value()*AssemMatX(i,j);
        for( int k=0; k<Sys->Ay.outerSize(); ++k)
          for (SpMatrix::InnerIterator it(Sys->Ay,k); it; ++it ) 
            Atotal.coeffRef(offset+it.row()*nBasis+i,offset+it.col()*nBasis+j) += it.value()*AssemMatY(i,j);
        for( int k=0; k<Sys->P.outerSize(); ++k)
          for (SpMatrix::InnerIterator it(Sys->P,k); it; ++it ) 
            Atotal.coeffRef(offset+it.row()*nBasis+i,offset+it.col()*nBasis+j) += it.value()*AssemMass(i,j);
            
        MatrixXd pos0 = pos.col(j);
        MatrixXd f = Sys->Force(pos0);
        for( int k=0; k<Sys->nEqn; k++)
          rhs(offset+k*nBasis+i) += AssemMass(i,j)*f(k);
      };  
    };    
    
    VectorXi ID(type);
    Mesh->tag_get_data(Mesh->ID,&vertices[0],type,&ID(0));
    
    int nBasisE = order+1;
    MatrixXi idx(type,nBasisE),idxN(type,nBasisE);
    if( (order == 2) && (type == 3) ) {
      idx << 0,1,3, 1,2,4, 2,0,5;
      idxN << 0,2,5, 1,0,3, 2,1,4;    
    };
    if( (order == 1) && (type == 3) ) {
      idx << 0,1, 1,2, 2,0;
      idxN << 0,2, 1,0, 2,1;    
    };
    if( (order == 2) && (type == 4) ) {
      idx <<  0,1,4, 1,2,5, 2,3,6, 3,0,7;
      idxN << 0,3,7, 1,0,4, 2,1,5, 3,2,6;    
    };
    if( (order == 1) && (type == 4) ) {
      idx <<  0,1, 1,2, 2,3, 3,0;
      idxN << 0,3, 1,0, 2,1, 3,2;    
    };
    
    // edge integrals
    for( int nn = 0; nn < type; nn++ ) {
      MatrixXd posE(3,nBasisE);
      for( int ix=0; ix < nBasisE; ix++ ) posE.col(ix) = pos.col(idx(nn,ix));
      
      MatrixXd normal = cross(posE.col(0)-posE.col(1));
      MatrixXd Am = Sys->Aminus(normal);
      MatrixXd AssemMassE = AssembleMassE(posE,order);
      
      Range neighbor;
      Mesh->get_adjacencies(&(vertices[nn]), 2, 2, true, neighbor );
      neighbor -= Range(*elem,*elem);
      
      int offsetN =  Mesh->get_ID(neighbor.front())*Sys->nEqn*nBasis;
            
      // interior edge
      if(neighbor.size() == 1 ) {
        const EntityHandle *connectN;
        int num_connectN;	
        Mesh->get_connectivity(neighbor.front(), connectN, num_connectN, true);
        VectorXi IDN(type);                                          
        Mesh->tag_get_data(Mesh->ID,connectN,num_connectN,&IDN(0));
        int pp = 0; while( pp < type && ID[nn] != IDN[pp]) { pp++; }
        
        // line integrals
        for( int i=0; i<nBasisE; i++ )
          for( int j=0; j<nBasisE; j++ )
            for( int k=0; k<Sys->nEqn; k++ )
              for( int l=0; l<Sys->nEqn; l++ ) {
                Atotal.coeffRef(offset+k*nBasis+idx(nn,i),offset+l*nBasis+idx(nn,j)) += 0.5*Am(k,l)*AssemMassE(i,j);
                Atotal.coeffRef(offset+k*nBasis+idx(nn,i),offsetN+l*nBasis+idxN(pp,j)) -= 0.5*Am(k,l)*AssemMassE(i,j);
              };
      };
      
      // boundary edge
      if(neighbor.size() == 0 ) {
        Range edge;
        Mesh->get_adjacencies(&(vertices[nn]), 2, 1, true, edge );
        int boundaryflag = Mesh->tag_get_integer(Mesh->BoundaryID,edge.front());
      
        MatrixXd normals(3,3);
        Mesh->tag_get_data(Mesh->BoundaryNormal,&(vertices[nn]),2,&normals(0)); 
        normals.col(2) = normal;
        
        // line integrals
        for( int i=0; i<nBasisE; i++ ) {
          for( int j=0; j<nBasisE; j++ ) { 
            normal = NormalApprox(normals,j,nBasisE);
            MatrixXd pos0 = posE.col(j);
            Sys->setBCData(pos0,normal,boundaryflag);
            
            MatrixXd AmBC = 0.5*Am*Sys->XBC;
            for( int k=0; k<Sys->nEqn; k++ )
              for( int l=0; l<Sys->nEqn; l++ )
                Atotal.coeffRef(offset+k*nBasis+idx(nn,i),offset+l*nBasis+idx(nn,j)) += AmBC(k,l)*AssemMassE(i,j);	     
              
            MatrixXd AmBCrhs = 0.5*Am*Sys->BCrhs;	    
            for( int k=0; k<Sys->nEqn; k++)
              rhs(offset+k*nBasis+idx(nn,i)) += -AssemMassE(i,j)*AmBCrhs(k);	    
          };  
        };
      };      
    };
    
  };
  stopTimer();
  
  startTimer(" compress sparse matrix\n"); 
  Atotal.prune(1,1e-13);  
  Atotal.makeCompressed();  
  stopTimer();
  
  nnz = Atotal.nonZeros();
  cout << "Non zeros: " << nnz << endl;
  cout << "Non zeros per row: " << 1.0*nnz/nTot << endl;
  
};

// normal approximation at boundary edge Lagrange points
MatrixXd Numerics::NormalApprox( const MatrixXd& normals, int j, int nBasisE )
{
  MatrixXd normal;
  double diff = (normals.col(0)-normals.col(1)).norm();
  
  if( diff < 0.6 ) {
    if( j == 2 ) return( (0.5*normals.col(0)+0.5*normals.col(1)).normalized() );
    else return( normals.col(j) );
  } else {
    return( normals.col(2) );
  };
};

// generate 'cross-vector' in 2D 
MatrixXd Numerics::cross( const MatrixXd& vector )
{
  MatrixXd cross(3,1);
  cross << -vector(1,0),vector(0,0),0;
  return( cross.normalized() );
};

// linear solve using direct sparse solver with LU decomposition
void Numerics::LinearSolve()
{
  INT iparm[64];
  double dparm[64];
  void *pt[64];
  for(int i=0; i<64; i++ ) pt[i]=0;

  startTimer(">> LU decomposition\n"); 

  PardisoSolve(11, Atotal.rows(), (INT *)Atotal.outerIndexPtr(), (INT *)Atotal.innerIndexPtr(), Atotal.valuePtr(), 
               Atotal.nonZeros(), &rhs(0), &sol(0), pt, -11, iparm, dparm);
  stopTimer(); 
  startTimer("");
  PardisoSolve(11, Atotal.rows(), (INT *)Atotal.outerIndexPtr(), (INT *)Atotal.innerIndexPtr(), Atotal.valuePtr(), 
               Atotal.nonZeros(), &rhs(0), &sol(0), pt, 11, iparm, dparm);
  stopTimer(); 
  startTimer("");
  PardisoSolve(11, Atotal.rows(), (INT *)Atotal.outerIndexPtr(), (INT *)Atotal.innerIndexPtr(), Atotal.valuePtr(), 
               Atotal.nonZeros(), &rhs(0), &sol(0), pt, 22, iparm, dparm);  
  stopTimer(); 

  startTimer(" Solve\n");  
  PardisoSolve(11, Atotal.rows(), (INT *)Atotal.outerIndexPtr(), (INT *)Atotal.innerIndexPtr(), Atotal.valuePtr(), 
               Atotal.nonZeros(), &rhs(0), &sol(0), pt, 33, iparm, dparm);    
  nIter = iparm[6];
  stopTimer(); 

  // free memory
  PardisoSolve(11, Atotal.rows(), (INT *)Atotal.outerIndexPtr(), (INT *)Atotal.innerIndexPtr(), Atotal.valuePtr(), 
               Atotal.nonZeros(), &rhs(0), &sol(0), pt, -1, iparm, dparm);    

  startTimer(" Compute residuum\n");  
  rhs = Atotal*sol-rhs;
  res = rhs.norm();
  stopTimer();
  
  cout << "Residuum: " << res << endl;
};

// writes/attaches solution to mesh object (element center values and vertex values)
void Numerics::AttachSolution()
{
  int length = Sys->nEqn;
  if( length > 9 ) length = 9; // restrict to essential fields

  Mesh->tag_get_handle("Solution Vertices", length, MB_TYPE_DOUBLE, SolutionVerts, MB_TAG_DENSE | MB_TAG_CREAT );
  for (Range::iterator vert = Mesh->Vertices.begin(); vert != Mesh->Vertices.end(); vert++) {
    int idx = Mesh->get_ID(*vert);
    Range Elems;
    MatrixXd solVec;
    solVec = MatrixXd::Zero(length,1);
    
    Mesh->get_adjacencies(&(*vert), 1, 2, true, Elems );
    for (Range::iterator elem = Elems.begin(); elem != Elems.end(); elem++) {
      const EntityHandle *vertices;
      int type;
      Mesh->get_connectivity(*elem, vertices, type, true);
      VectorXi ID(type);
      Mesh->tag_get_data(Mesh->ID,vertices,type,&ID(0));
      int pp = 0; while( pp < 3 && ID[pp] != idx) { pp++; }
      int offset = Mesh->get_ID(*elem)*Sys->nEqn*nBasis;
      
      for( int i=0; i<length; i++ ) 
        solVec(i) += sol(offset+i*nBasis+pp) / Elems.size();       
    }; 
    Mesh->tag_set_data(SolutionVerts,&(*vert),1,&solVec(0));      
  };  
  
  Mesh->tag_get_handle("Solution Elements", length, MB_TYPE_DOUBLE, SolutionElems, MB_TAG_DENSE | MB_TAG_CREAT );  
  for (Range::iterator elem = Mesh->Elements.begin(); elem != Mesh->Elements.end(); elem++) { 
    int offset =  Mesh->get_ID(*elem)*Sys->nEqn*nBasis;
    const EntityHandle *vertices;
    int type;
    Mesh->get_connectivity(*elem, vertices, type, true);
    MatrixXd solVec;
    solVec = MatrixXd::Zero(length,1);
    
    MatrixXd weights(nBasis,1);
    if( (order == 2) && (type == 3) ) weights << -1.0/9, -1.0/9, -1.0/9, 4.0/9, 4.0/9, 4.0/9;
    if( (order == 1) && (type == 3) ) weights << 1.0/3, 1.0/3, 1.0/3;
    if( (order == 2) && (type == 4) ) weights << 1.0/36, 1.0/36, 1.0/36, 1.0/36, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 4.0/9;
    if( (order == 1) && (type == 4) ) weights << 1.0/4, 1.0/4, 1.0/4, 1.0/4;
    for( int i=0; i<length; i++ )
      for( int j=0; j<nBasis; j++ ) 
        solVec(i) += weights(j)*sol(offset+i*nBasis+j);
    Mesh->tag_set_data(SolutionElems,&(*elem),1,&solVec(0));      
  };  
};

