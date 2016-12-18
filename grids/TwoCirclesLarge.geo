/*************************************************************************\
 ComplexChannel.geo  - geometry script for two con-centric circles.

This is to be used with the Gmsh program
http://gmsh.info/
\*************************************************************************/

// This value gives the global element size factor (lower -> finer mesh)
Mesh.CharacteristicLengthFactor = 0.5;

// This fixes the coloring to element type (better visibility)
Mesh.ColorCarousel = 0;

// Points are defined by {x,y,z,fac}, where fac gives a 
// relativ local refinement factor.

// three points to define two circles
Point(1) = {0, 0, 0, 1};
Point(2) = {0.5, 0, 0, 1};
Point(3) = {2, 0, 0, 1};

// Circles (or rather arcs) are constructed from a start point, 
// a center point and an end point. 

// two con-centric circles
Circle(1) = {2, 1, 2};
Circle(2) = {3, 1, 3};

// Polygons and more general line loops are constructed from simple 
// lines. A minus sign substracts the object.

Line Loop(3) = {1};
Line Loop(4) = {2};

// Surfaces are build similarly, second entries are considered holes.

Plane Surface(5) = {4,3};

// Finally, we assign labels (10,20) to the lines or group 
// of lines. Gmsh seems to require a label for the surface, too.

Physical Line(10) = {1};
Physical Line(20) = {2};
Physical Surface(100) = {5};
