// ###############################
// ###############################
// 2D cantilever beam
// Authors: Serafeim Bakalakos
// ###############################
// ###############################


// *******************************
// User defined parameters
// *******************************

// Parametric dimensions
minX = 0.0;
minY = 0.0;
maxX = 4.0;
maxY = 20.0;

// Mesh size
lc =  1.0;

// Element type: choose between tri3, quad4, tri6, quad8, quad9 by entering the appropriate number of nodes
elementNodes = 8;

// *******************************
// Script body
// *******************************

// Boundary points
Point(1) =  {minX, minY, 0.0, lc};
Point(2) =  {maxX, minY, 0.0, lc};
Point(3) =  {maxX, maxY, 0.0, lc};
Point(4) =  {minX, maxY, 0.0, lc};

// Boundary lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// Create individual surfaces from lines
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};


// Mesh options (I cannot get ElseIf to work lol)
If(elementNodes == 3)
	Mesh.ElementOrder = 1;
	Mesh.SecondOrderIncomplete = 0;
EndIf
If(elementNodes == 4)
	Mesh.ElementOrder = 1;
	Mesh.SecondOrderIncomplete = 0;
	Recombine Surface{1};
EndIf
If(elementNodes == 6)
	Mesh.ElementOrder = 2;
	Mesh.SecondOrderIncomplete = 0;
	Mesh.ElementOrder = 2;
EndIf
If(elementNodes == 8)
	Mesh.ElementOrder = 2;
	Mesh.SecondOrderIncomplete = 1;
	Recombine Surface{1};
EndIf
If(elementNodes == 9)
	Mesh.ElementOrder = 2;
	Mesh.SecondOrderIncomplete = 0;
	Recombine Surface{1};
EndIf


// Physical entities (only elements belonging to them will be saved)
Physical Surface(1) = {1};
