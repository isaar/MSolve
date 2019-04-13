// 
// Geometry file for structural member with fillet
// 
// Written by: Manolis Georgioudakis <geoem@mail.ntua.gr>


// Characteristic length
lc =  10.0;

// Parametric Dimensions b1, b2, r3
b1 = 75.0;
b2 = 75.0;
r3 = 20.0;

// Auxiliary parameters
bottom = 375.0;
h1 = (bottom - 2*r3 -b2)/2.0;
b3 = 75.0;

// Points definitions
Point(1) =  {           0.0,      0.0, 0.0, 2*lc};
Point(2) =  {           0.0,       b1, 0.0, 2*lc};
Point(3) =  {       h1-3*lc,       b1, 0.0, lc};
Point(4) =  {            h1,       b1, 0.0, lc/4};
Point(5) =  {            h1,    b1+r3, 0.0, lc};
Point(6) =  {         h1+r3,    b1+r3, 0.0, lc/4};
Point(7) =  {         h1+r3,    b1+b3, 0.0, lc/2};
Point(8) =  {      h1+r3+b2,    b1+b3, 0.0, lc/2};
Point(9) =  {      h1+r3+b2,    b1+r3, 0.0, lc/4};
Point(10) = {      bottom-h1,    b1+r3, 0.0, lc};
Point(11) = {      bottom-h1,       b1, 0.0, lc/4};
Point(12) = { bottom-h1+3*lc,       b1, 0.0, lc};
Point(13) = {         bottom,       b1, 0.0, 2*lc};
Point(14) = {         bottom,      0.0, 0.0, 2*lc};
Point(15) = { bottom-h1+3*lc,        0, 0.0, lc};
Point(16) = {       h1-3*lc,        0, 0.0, lc};

// Create boundary
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Circle(4) = {4,5,6};
Line(5) = {6,7};
Line(6) = {7,8};
Line(7) = {8,9};
Circle(8) = {9,10,11};
Line(9) = {11,12};
Line(10) = {12,13};
Line(11) = {13,14};
Line(12) = {14,15};
Line(13) = {15,16};
Line(14) = {16,1};
Line(15) = {16,3};
Line(16) = {6,9};
Line(17) = {12,15};

// Create individual surfaces from lines
Line Loop(1) = {1,2,-15,14};
Line Loop(2) = {15,3,4,16,8,9,17,13};
Line Loop(3) = {-17,10,11,12};
Line Loop(4) = {5,6,7,-16};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

// Recombine the triangles into quads:
Recombine Surface{1,2,3,4};

// Physical entities (for which elements will be saved)
Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};