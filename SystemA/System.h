/*************************************************************************\
 System.h  - declaration of system class. 
\*************************************************************************/
  
struct solution;

class System 
{
  public:
    System( int argc, const char **argv );
    int nEqn, nBC;
    SpMatrix Ax, Ay, P, BC;
    VectorXi OddVar, TensorDegree;
    MatrixXd XBC, BCrhs;
    void setBCData( MatrixXd& pos, MatrixXd& normal, int boundaryID );
    void setEQNData();
    
    MatrixXd Force( MatrixXd& pos );
    
    MatrixXd Aminus( const MatrixXd& normal );
    
    double A0, A1, A2;
    double tau, chi, zeta, kappa, theta0, theta1, uW, v0, eps;
    string SystemText;
    
    string PostProcess( MeshCore *Mesh );
    string InfoText();
    
  private:
    MatrixXd X0, Aminus1D;
    int outputFlag;
    
    string FileOutput( MeshCore *Mesh );
    string ComputeError( MeshCore *Mesh );
    SpMatrix Projector( const MatrixXd& normal );
    SpMatrix invProjector( const MatrixXd& normal );
    
    solution *solX;
    void setupExact();    
};

