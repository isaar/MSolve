// Parameters
L = 5.0;
h = 0.5;
lc = 0.1;

// Points
Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {L, 0.0, 0.0, lc};
Point(3) = {L, h, 0.0, lc};
Point(4) = {0.0, h, 0.0, lc};

// Boundary
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// Surface
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};