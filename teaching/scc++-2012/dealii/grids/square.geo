nx = 10;
ny = 10;

// Corner points of unit square
Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {1, 1, 0};
Point(4) = {0, 1, 0};

// Sides of unit square
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Transfinite Surface(1) = {1,2,3,4};
Recombine Surface(1);

// Specify number of points on each line
Transfinite Line{1,3} = nx;
Transfinite Line{2,4} = ny;

// Here we give a tag = 0 for each line
Physical Line(0) = {1,2,3,4};
Physical Surface(0) = {1};
