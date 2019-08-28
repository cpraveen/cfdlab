lc = 0.02;
xmin = -1.0;
xmax = +1.0;
ymin = -1.0;
ymax = +1.0;

Point(1) = {xmin,ymin,0,lc};
Point(2) = {xmax,ymin,0,lc};
Point(3) = {xmax,ymax,0,lc};
Point(4) = {xmin,ymax,0,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Physical Surface(100) = {1};
Physical Line(200) = {1};
Physical Line(300) = {2};
Physical Line(400) = {3};
Physical Line(500) = {4};
