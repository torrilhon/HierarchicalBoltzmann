/*************************************************************************\
 ComplexChannel.geo  - geometry script for a curved channel with 
                       three obstacles.
This is to be used with the Gmsh program
http://gmsh.info/
\*************************************************************************/

// This value gives the global element size factor (lower -> finer mesh)
Mesh.CharacteristicLengthFactor = 0.1;

// This fixes the coloring to element type (better visibility)
Mesh.ColorCarousel = 0;

// Points are defined by {x,y,z,fac}, where fac gives a 
// relativ local refinement factor.

// control points for the channel
Point(1) = {-2, 1, 0, 1};
Point(2) = {-2, 0, 0, 1};
Point(3) = {0, 0, 0, 1};
Point(4) = {0, -2, 0, 1};
Point(5) = {1, -2, 0, 1};
Point(6) = {1, 1, 0, 1};

//control points obstacle 1
Point(11) = {-1.4, 0.35, 0, 0.3};
Point(12) = {-1, 0.5, 0, 0.3};
Point(13) = {-0.6, 0.35, 0, 0.3};
Point(14) = {-0.6, 0.6, 0, 0.3};
Point(15) = {-1, 0.75, 0, 0.3};
Point(16) = {-1.4, 0.75, 0, 0.3};

//control points obstacle 2
Point(21) = {-0.1, -0.2, 0, 0.3};
Point(22) = {0.25, -0.25, 0, 0.3};
Point(23) = {0.4, -0.1, 0, 0.3};
Point(24) = {0.1, 0, 0, 0.3};
Point(25) = {-0.1, 0.3, 0, 0.3};
Point(26) = {-0.4, 0.1, 0, 0.3};

//control points obstacle 3
Point(31) = {0.55, -0.8, 0, 0.3};
Point(32) = {0.2, -1.2, 0, 0.3};
Point(33) = {0.6, -1.3, 0, 0.3};
Point(34) = {0.8, -0.9, 0, 0.3};
Point(35) = {0.7, -0.4, 0, 0.3};
Point(36) = {0.5, -0.4, 0, 0.3};

// Lines are constructed from the points. Bsplines take control 
// points. Careful: lines follow a direction.

// channel shape
Line(1) = {1, 2};
BSpline(2) = {2, 3, 4};
Line(3) = {4, 5};
BSpline(4) = {5, 6, 1};

// obstacle 1
BSpline(5) = {11, 12, 13, 14, 15, 16, 11};

// obstacle 2
BSpline(6) = {21, 22, 23, 24, 25, 26, 21};

// obstacle 3
BSpline(7) = {31, 32, 33, 34, 35, 36, 31};

// Polygons and more general line loops are constructed from simple 
// lines. A minus sign substracts the object.

Line Loop(8) = {1, 2, 3, 4};
Line Loop(9) = {5};
Line Loop(10) = {6};
Line Loop(11) = {7};

// Surfaces are build similarly, second entries are considered holes.

Plane Surface(12) = {8, 9, 10, 11};

// Finally, we assign labels (10,20,30,40) to the lines or group 
// of lines. Gmsh seems to require a label for the surface, too.

Physical Line(10) = {3};    // inflow
Physical Line(20) = {1};    // outflow 
Physical Line(30) = {2, 4}; // wall
Physical Line(40) = {5};    // obstacle
Physical Line(50) = {6};    // obstacle
Physical Line(60) = {7};    // obstacle
Physical Surface(100) = {12};
