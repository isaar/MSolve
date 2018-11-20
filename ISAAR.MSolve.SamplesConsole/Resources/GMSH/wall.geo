// Mesh sizes
meshSizeExternal = 0.20;
meshSizeInternal = 0.15;

// Wall params
wallLength = 6.0;
wallHeight = 3.5;

// Door params
doorCenterX = wallLength / 2;
doorWidth = 0.9144;
doorHeight = 2.032;

// Window 1 params
window1CenterX = wallLength / 5;
window1TopY = doorHeight;
window1Width = 0.6096;
window1Height = 0.9144;

// Window 2 params
window2CenterX = wallLength - window1CenterX;
window2TopY = doorHeight;
window2Width = window1Width;
window2Height = window1Height;

// Auxiliary quantities
doorLeftX = doorCenterX - doorWidth / 2;
doorRightX = doorCenterX + doorWidth / 2;
window1MinX = window1CenterX - window1Width / 2;
window1MinY = window1TopY - window1Height;
window1MaxX = window1CenterX + window1Width / 2;
window1MaxY = window1TopY;
window2MinX = window2CenterX - window2Width / 2;
window2MinY = window2TopY - window2Height;
window2MaxX = window2CenterX + window2Width / 2;
window2MaxY = window2TopY;

// Points definitions
// External
Point(1) = { 0.0, 0.0, 0.0, meshSizeExternal };
Point(2) = { doorLeftX, 0.0, 0.0, meshSizeExternal };
Point(3) = { doorLeftX, doorHeight, 0.0, meshSizeExternal };
Point(4) = { doorRightX, doorHeight, 0.0, meshSizeExternal };
Point(5) = { doorRightX, 0.0, 0.0, meshSizeExternal };
Point(6) = { wallLength, 0.0, 0.0, meshSizeExternal };
Point(7) = { wallLength, wallHeight, 0.0, meshSizeExternal };
Point(8) = { 0.0, wallHeight, 0.0, meshSizeExternal };
// Internal
Point(9) =  { window1MinX, window1MinY, 0.0, meshSizeInternal };
Point(10) = { window1MaxX, window1MinY, 0.0, meshSizeInternal };
Point(11) = { window1MaxX, window1MaxY, 0.0, meshSizeInternal };
Point(12) = { window1MinX, window1MaxY, 0.0, meshSizeInternal };
Point(13) = { window2MinX, window2MinY, 0.0, meshSizeInternal };
Point(14) = { window2MaxX, window2MinY, 0.0, meshSizeInternal };
Point(15) = { window2MaxX, window2MaxY, 0.0, meshSizeInternal };
Point(16) = { window2MinX, window2MaxY, 0.0, meshSizeInternal };

// Create boundary
//External
Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 3, 4 };
Line(4) = { 4, 5 };
Line(5) = { 5, 6 };
Line(6) = { 6, 7 };
Line(7) = { 7, 8 };
Line(8) = { 8, 1 };
//Internal 1
Line(9)  = {  9, 10 };
Line(10) = { 10, 11 };
Line(11) = { 11, 12 };
Line(12) = { 12,  9 };
//Internal 2
Line(13) = { 13, 14 };
Line(14) = { 14, 15 };
Line(15) = { 15, 16 };
Line(16) = { 16, 13 };

// Create individual surfaces from lines
Line Loop(1) = { 1, 2, 3, 4, 5, 6, 7, 8 };
Line Loop(2) = { 9, 10, 11, 12} ;
Line Loop(3) = { 13, 14, 15, 16 };
Plane Surface(1) = { 1, 2, 3 };

// Recombine the triangles into quads:
Recombine Surface{1};

// Physical entities (for which elements will be saved)
Physical Surface(1) = { 1 };