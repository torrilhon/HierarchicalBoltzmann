/*************************************************************************\
 Mesh.h  - declares derived mesh class with some specific functions.
           Definitions are in Mesh.cpp.
 This is based on the MOAB library, see 
 http://sigma.mcs.anl.gov/moab-library/
\*************************************************************************/
#include "moab/Core.hpp"
#include "moab/Range.hpp"
using namespace moab;

class MeshCore : public Core  
{
  public:
    MeshCore();
    MeshCore( int argc, const char **argv );
    void setup();
    int get_ID( const EntityHandle& entity_handle );
    int tag_get_integer( const Tag tag_handle, const EntityHandle& entity_handle );
    string InfoText();
    
    double hMax;
    int Type;
    Range Vertices, Edges, Elements, BCVertices;
    Tag ID, BoundaryTag, BoundaryNormal, BoundaryID;
    string MeshFile;
    
  private:
    MatrixXd NormalApprox( const MatrixXd& P0, const MatrixXd& Pts );
};
