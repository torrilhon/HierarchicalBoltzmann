/*************************************************************************\
 Mesh.cpp  - reads and initializes the mesh
 This is based on the MOAB library, see 
 http://sigma.mcs.anl.gov/moab-library/
\*************************************************************************/
#include <stdio.h>
#include <iostream>
#include <math.h>
using namespace std;

#include "EigenSetup.h"
#include "Mesh.h"
#include "Tools.h"


// writes logfile message
string MeshCore::InfoText()
{
  return(strprintf("MeshFile = '%s', type = %d, hMax = %.6f, nElems = %d, nVerts = %d", 
                    MeshFile.c_str(), Type, hMax, Elements.size(), Vertices.size() ) );
};

// default constructor
MeshCore::MeshCore() : Core()
{
};

// constructor with command line arguments
MeshCore::MeshCore( int argc, const char **argv ) : Core()
{
  MeshFile = argv[1];

  setup();
};

// reads and initializes the mesh
void MeshCore::setup() 
{
  int bcFlagCount = 0;
    
  startTimer(">> Reading Mesh\n"); {
    cout << MeshFile << "\n";   
    int ret = load_mesh(MeshFile.c_str());
    if( ret != 0 ) {
      cout << " Mesh file not found !!!!!\n";
      exit(0);
    };
  } stopTimer();
  
  
  startTimer(" generate edges \n"); {
    get_entities_by_type(0,MBVERTEX,Vertices);
    
    for (Range::iterator vert = Vertices.begin(); vert != Vertices.end(); vert++) {
      Range adjs;
      get_adjacencies(&(*vert), 1, 1, true, adjs);
      if( adjs.size() == 0 ) delete_entities(&(*vert),1);
    };
    get_entities_by_type(0,MBEDGE,Edges);
    Range Triangles, Quads;
    get_entities_by_type(0,MBTRI,Triangles);
    get_entities_by_type(0,MBQUAD,Quads);
    Elements.insert(Triangles.begin(),Triangles.end());
    Elements.insert(Quads.begin(),Quads.end());
    if( Quads.size() > 0 ) Type = 4;
    if( Triangles.size() > 0 ) Type = 3;
  } stopTimer();
  
  startTimer(" compute and store boundary normal, edge length, id's, etc on the mesh \n"); {
    // we delete most of the existing tags (note: GLOBAL_ID seems needed for 'list_entities()')    
    Tag Tmp;
    if(tag_get_handle("GLOBAL_ID",Tmp)==0) tag_delete(Tmp);
    if(tag_get_handle("NEUMANN_SET",Tmp)==0) tag_delete(Tmp);
    if(tag_get_handle("DIRICHLET_SET",Tmp)==0) tag_delete(Tmp);
    if(tag_get_handle("GEOM_DIMENSION",Tmp)==0) tag_delete(Tmp);

    // MATERIAL_SET flags the physical domains (from gmsh), we distribute this to the edges
    // DomainFlag: flag of the boundary edges
    tag_get_handle("MATERIAL_SET",Tmp);
    tag_get_handle("Boundary ID", 1, MB_TYPE_INTEGER, BoundaryID, MB_TAG_SPARSE | MB_TAG_CREAT ); 

    Range Sets;
    get_entities_by_type_and_tag(0, MBENTITYSET, &Tmp, NULL, 1, Sets, Interface::UNION);
    for( Range::iterator set = Sets.begin(); set != Sets.end(); set++)  {
      int value;
      Range SetEdges;
      tag_get_data(Tmp, &(*set), 1, &value);
      get_entities_by_dimension(*set, 1, SetEdges);
      for (Range::iterator edge = SetEdges.begin(); edge != SetEdges.end(); edge++) { 
        tag_set_data(BoundaryID,&(*edge),1,&value);
        bcFlagCount++;
      };
    };

    // we create our own ID tag 
    // ID: to number vertices with SparseMatrix indices
    tag_get_handle("Entity ID", 1, MB_TYPE_INTEGER, ID, MB_TAG_DENSE | MB_TAG_CREAT );   
    get_entities_by_type(0,MBVERTEX,Vertices);
    int count = 0, zero = 0, one[3] = {1,1,1};
    for (Range::iterator vert = Vertices.begin(); vert != Vertices.end(); vert++) {
      tag_set_data(ID,&(*vert),1,&count);
      count++;
    };
    count = 0;
    for (Range::iterator elem = Elements.begin(); elem != Elements.end(); elem++) {
      tag_set_data(ID,&(*elem),1,&count);
      count++;
    };
    
    // BoundaryTag: to identify boundary entities
    tag_get_handle("Boundary Flag", 1, MB_TYPE_INTEGER, BoundaryTag, MB_TAG_DENSE | MB_TAG_CREAT ); 
    for (Range::iterator edge = Edges.begin(); edge != Edges.end(); edge++) { 
      Range elements, vertices;
      get_adjacencies(&(*edge), 1, 2, true, elements);
      if( elements.size() == 1 ) {
        tag_set_data(BoundaryTag,&(*edge),1,&one);  
        tag_set_data(BoundaryTag,elements,&one);
        get_adjacencies(&(*edge), 1, 0, true, vertices);    
        tag_set_data(BoundaryTag,vertices,&one);  
      };
    };
    
    // BoundaryNormal: to accommodate the boundary normal
    tag_get_handle("Boundary Normal", 3, MB_TYPE_DOUBLE, BoundaryNormal, MB_TAG_SPARSE | MB_TAG_CREAT );
    void *onePtr = one;
    get_entities_by_type_and_tag(0,MBVERTEX,&BoundaryTag,&onePtr,1,BCVertices);   
    for (Range::iterator vert = BCVertices.begin(); vert != BCVertices.end(); vert++) {
      MatrixXd pos0(3,1), normal(3,1);
      Range elements, vertices, vertsOUT, vertsIN;
      get_coords(&(*vert),1,&pos0(0));
      get_adjacencies(&(*vert), 1, 2, true, elements);
      
      if( elements.size() > 1 ) // vertex is true point...
        get_connectivity(elements, vertices, true); // only corners 
      else // vertex is edge mid-point...
        get_connectivity(elements, vertices, false); // not only corners 
        
      for (Range::iterator it = vertices.begin(); it != vertices.end(); it++) 
        if(tag_get_integer(BoundaryTag,*it)) vertsOUT.insert(*it);  
        else vertsIN.insert(*it);       
          
      MatrixXd posOUT(3,vertsOUT.size()); 
      MatrixXd posIN(3,vertsIN.size());
      get_coords(vertsOUT,&posOUT(0));
      get_coords(vertsIN,&posIN(0));
      normal = NormalApprox(pos0,posOUT);
      double outside = -(normal.transpose()*(posIN.rowwise().mean()-pos0.col(0)))(0);
      if( outside < 0 ) normal = -normal;  
      tag_set_data(BoundaryNormal,&(*vert),1,&normal(0));      
    };
    
    // compute maximal edge length hMax
    hMax = 0.0;
    for (Range::iterator edge = Edges.begin(); edge != Edges.end(); edge++) { 
      Range elements, vertices;
      MatrixXd pos(3,2);
      get_connectivity(&(*edge), 1, vertices, true);
      get_coords(vertices,&pos(0));
      double length = (pos.col(1)-pos.col(0)).norm();
      if( length > hMax ) hMax = length;
    };
  } stopTimer();
  
  cout << "Type of elements " << Type << "\n";
  cout << "Number of corner vertices " << Vertices.size() << "\n";
  cout << "Number of BC vertices " << BCVertices.size() << "\n";
  cout << "Number of edges " << Edges.size() << "\n";
  cout << "Number of flagged BC edges " << bcFlagCount << "\n";
  cout << "Maximal edge length " << hMax << "\n";
  cout << "Number of elements " << Elements.size() << "\n";
    
};

// returns value of ID-tag at at entity 'entity_handle' directly
int MeshCore::get_ID( const EntityHandle& entity_handle )
{
  int ret = -1;
  tag_get_data(ID,&entity_handle,1,&ret);
  return( ret );
};

// returns value of an integer tag with 'tag_handle' at entity 'entity_handle'
int MeshCore::tag_get_integer( const Tag tag_handle, const EntityHandle& entity_handle )
{
  int ret;
  tag_get_data(tag_handle,&entity_handle,1,&ret);
  return( ret );
};


// Approximation of a normal vector from three points
// input: position of normal P0, three points Pts(3,3) defining a circle
// the normal is computed from the direction between the midpoint of
// the circle and P0, a line as degenerated circle is included.
MatrixXd MeshCore::NormalApprox( const MatrixXd& P0, const MatrixXd& Pts )
{
  MatrixXd P1 = Pts.block(0,0,3,1), P2 = Pts.block(0,1,3,1), P3 = Pts.block(0,2,3,1); 
  double a1 = P1(0), b1 = P1(1), a2 = P2(0), b2 = P2(1), a3 = P3(0), b3 = P3(1);
  double norm1 = P1.squaredNorm(), norm2 = P2.squaredNorm(), norm3 = P3.squaredNorm();
  double det = a3*(b2 - b1) + a2*(b1 - b3) + a1*(-b2 + b3);
  MatrixXd P(3,1), normal(3,1);
  
  if( fabs(det) > 1e-11 ) {
    P(0) = 0.5*(b3*(norm1-norm2) + b1*(norm2-norm3) + b2*(norm3-norm1))/det;
    P(1) = 0.5*(a3*(norm2-norm1) + a2*(norm1-norm3) + a1*(norm3-norm2))/det;
    P(2) = 0.0;
    normal = (P0 - P).normalized();
  } else {
    P = (P3 - P1).normalized();
    normal << -P(1), P(0), 0.0;
  };
  
  return(normal);    
};




