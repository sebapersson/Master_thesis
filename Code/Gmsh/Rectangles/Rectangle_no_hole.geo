lc = 0.1;

// Outer rectangle points 
Point(1) = {-3.0, -3.0, 0, lc};
Point(2) = {3.0, -3.0, 0, lc};
Point(3) = {3.0, 3.0, 0, lc};
Point(4) = {-3.0, 3.0, 0, lc};

// Lines between the points in the rectangle
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Loop over the lines 
Line Loop(5) = {1, 2, 3, 4};

// Add a plane surface for the rectangle
Plane Surface(1) = {5};

// Make the boundaries to physical lines
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};

// Make the surface to a physical surface
Physical Surface(1) = {1};
