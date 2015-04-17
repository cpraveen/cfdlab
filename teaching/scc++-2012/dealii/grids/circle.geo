Mesh.RecombineAll=1;           //recombine all defined surfaces
Mesh.Algorithm=8;              //delquad mesher
Mesh.RecombinationAlgorithm=1; //blossom

r = 1;
s = r/2;

n1 = 5;
n2 = 5;

Point(1) = { 0,  0, 0};
Point(2) = { s,  0, 0};
Point(3) = { 0,  s, 0};
Point(4) = {-s,  0, 0};
Point(5) = { 0, -s, 0};
Point(6) = { r,  0, 0};
Point(7) = { 0,  r, 0};
Point(8) = {-r,  0, 0};
Point(9) = { 0, -r, 0};

Line(1) = {2,3};
Line(2) = {3,4};
Line(3) = {4,5};
Line(4) = {5,2};
Circle(5) = {6,1,7};
Circle(6) = {7,1,8};
Circle(7) = {8,1,9};
Circle(8) = {9,1,6};
Line(9) = {2,6};
Line(10) = {3,7};
Line(11) = {4,8};
Line(12) = {5,9};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Transfinite Surface(1) = {2,3,4,5};

Line Loop(2) = {9,5,-10,-1};
Plane Surface(2) = {2};
Transfinite Surface(2) = {2,6,7,3};

Line Loop(3) = {10,6,-11,-2};
Plane Surface(3) = {3};
Transfinite Surface(3) = {3,7,8,4};

Line Loop(4) = {11,7,-12,-3};
Plane Surface(4) = {4};
Transfinite Surface(4) = {4,8,9,5};

Line Loop(5) = {12,8,-9,-4};
Plane Surface(5) = {5};
Transfinite Surface(5) = {5,9,6,2};

// The following are not required since we have
// specified global option to recombine all surfaces
//Recombine Surface(1);
//Recombine Surface(2);
//Recombine Surface(3);
//Recombine Surface(4);
//Recombine Surface(5);

Transfinite Line{1,2,3,4,5,6,7,8} = n1;
Transfinite Line{9,10,11,12} = n2;

Physical Line(0) = {5,6,7,8};
Physical Surface(0) = {1,2,3,4,5};
