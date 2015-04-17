Mesh.RecombineAll=1;           //recombine all defined surfaces
Mesh.Algorithm=8;              //delquad mesher
Mesh.RecombinationAlgorithm=1; //blossom

r = 1;

n = 50;

cl = (2*Pi*r)/n;

Point(1) = { 0,  0, 0, cl};
Point(2) = { r,  0, 0, cl};
Point(3) = { 0,  r, 0, cl};
Point(4) = {-r,  0, 0, cl};
Point(5) = { 0, -r, 0, cl};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Physical Surface(0) = {1};
Physical Line(0) = {1,2,3,4};
