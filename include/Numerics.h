/*************************************************************************\
 Numerics.h  - declaration of numerics class. Most of the definitions 
               are in Numerics.cpp, some are in ElementAssembly.cpp.
\*************************************************************************/

class Numerics  
{
  public:
    Numerics( MeshCore *Mesh, System *Sys, int argc, const char **argv );
    
    void Assemblation();
    void LinearSolve();
    void AttachSolution();
    string InfoText();
    
    MeshCore *Mesh;
    System *Sys;
    Tag SolutionVerts, SolutionElems; 
    int nBasis, order, nTot, nnz, nIter;
    double res;
    
  private:
    MatrixXd AssembleDiffX( const MatrixXd& pos, int order, int type );
    MatrixXd AssembleDiffY( const MatrixXd& pos, int order, int type );
    MatrixXd AssembleMass( const MatrixXd& pos, int order, int type );
    MatrixXd AssembleMassE( const MatrixXd& pos, int order );
    void LagrangePts( MatrixXd& Pts, int order, int type );
    MatrixXd cross( const MatrixXd& vector );
    MatrixXd NormalApprox( const MatrixXd& normals, int j, int nBasisE );
    
    SpMatrixL Atotal;
    MatrixXd sol, rhs;
};
