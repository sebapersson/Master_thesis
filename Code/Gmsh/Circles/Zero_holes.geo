lc = 0.04;

R = 2.5;
// Outer rectangle points 
Point(1) = {0, -R, 0, lc};
Point(2) = {R, 0, 0, lc};
Point(3) = {0, R, 0, lc};
Point(4) = {-R, 0 , 0, lc};
Point(1000) = {0, 0, 0, lc};

// Lines between the points in the rectangle
Circle(1) = {1, 1000, 2};
Circle(2) = {2, 1000, 3};
Circle(3) = {3, 1000, 4};
Circle(4) = {4, 1000, 1};
Line Loop(1) = {1, 2, 3, 4};

// Build the surface and add line 
Physical Line(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Surface(1) = {1};