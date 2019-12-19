lc = 0.08;

// Outer rectangle points 
limitRec = 2.2;
Point(1) = {-limitRec, -limitRec, 0, lc};
Point(2) = {limitRec, -limitRec, 0, lc};
Point(3) = {limitRec, limitRec, 0, lc};
Point(4) = {-limitRec, limitRec, 0, lc};

// Making the points of the circles 
r = 0.25;

// Points for the first circle
Point(5) = {0.0, 0.0, 0.0, lc};
Point(6) = {r, 0.0, 0.0, lc};
Point(7) = {-r, 0.0, 0.0, lc};

// Second circle 
Point(8) = {0.40, 0.40, 0.0, lc};
Point(9) = {0.40 + r, 0.40, 0.0, lc};
Point(10) = {0.40 - r, 0.40, 0.0, lc};

// Third circle 
Point(11) = {0.55, -0.15, 0.0, lc};
Point(12) = {0.55 + r, -0.15, 0.0, lc};
Point(13) = {0.55 - r, -0.15, 0.0, lc};

// Fourth circle
Point(14) = {-0.15, -0.50, 0.0, lc};
Point(15) = {-0.15 + r, -0.50, 0.0, lc};
Point(16) = {-0.15 - r, -0.50, 0.0, lc};

// Fifth circle
Point(17) = {-0.55, 0.1, 0, lc};
Point(18) = {-0.55 + r, 0.1, 0, lc};
Point(19) = {-0.55 - r, 0.1, 0, lc};

// Lines between the points in the rectangle
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Lines for the first circle
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 6};

// Second circle lines 
Circle(7) = {9, 8, 10};
Circle(8) = {10, 8, 9};

// Third circle line 
Circle(9) = {12, 11, 13};
Circle(10) = {13, 11, 12};

// Furth circle loop 
Circle(11) = {15, 14, 16};
Circle(12) = {16, 14, 15};

// Fifth circle 
Circle(13) = {18, 17, 19};
Circle(14) = {19, 17, 18};

// Loop over the lines to make them into one line
Line Loop(9) = {1, 2, 3, 4};
Line Loop(10) = {5, 6};
Line Loop(11) = {7, 8};
Line Loop(12) = {9, 10};
Line Loop(13) = {11, 12};
Line Loop(14) = {13, 14};

// Add a plane surface for the rectangle
Plane Surface(1) = {9, 10, 11, 12, 13, 14};
Plane Surface(2) = {10};
Plane Surface(3) = {11};
Plane Surface(4) = {12};
Plane Surface(5) = {13};
Plane Surface(6) = {14};

// Make the boundaries to physical lines
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Line(5) = {5, 6};
Physical Line(6) = {7, 8};
Physical Line(7) = {9, 10};
Physical Line(8) = {11, 12};
Physical Line(9) = {13, 14};

// Make the surface to a physical surface
Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};





