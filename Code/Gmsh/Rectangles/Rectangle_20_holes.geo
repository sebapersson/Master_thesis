lc = 0.1;

// Outer rectangle points 
Point(1) = {-3.0, -3.0, 0, lc};
Point(2) = {3.0, -3.0, 0, lc};
Point(3) = {3.0, 3.0, 0, lc};
Point(4) = {-3.0, 3.0, 0, lc};

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

// Sixth circle 
x6 = -0.15;
y6 = 0.50;
Point(20) = {x6, y6, 0, lc};
Point(21) = {x6 + r, y6, 0, lc};
Point(22) = {x6 - r, y6, 0, lc};

// 7:th circle 
x7 = 0.40;
y7 = -0.65;
Point(23) = {x7, y7, 0, lc};
Point(24) = {x7 + r, y7, 0, lc};
Point(25) = {x7 - r, y7, 0, lc};

// 8:th circle 
x8 = -0.70;
y8 = -0.40;
Point(26) = {x8, y8, 0, lc};
Point(27) = {x8 + r, y8, 0, lc};
Point(28) = {x8 - r, y8, 0, lc};

// 9:th circle 
x9 = -0.50;
y9 = -0.90;
Point(29) = {x9, y9, 0, lc};
Point(30) = {x9 + r, y9, 0, lc};
Point(31) = {x9 - r, y9, 0, lc};

// 10:th circle 
x10 = 0.02;
y10 = -1.0;
Point(32) = {x10, y10, 0, lc};
Point(33) = {x10 + r, y10, 0, lc};
Point(34) = {x10 - r, y10, 0, lc};

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
// Sixth circle 
Circle(15) = {21, 20, 22};
Circle(16) = {22, 20, 21};
// 7:th circle 
Circle(17) = {24, 23, 25};
Circle(18) = {25, 23, 24};
// 8:th circle 
Circle(19) = {28, 26, 27};
Circle(20) = {27, 26, 28};
// 9:th circle 
Circle(21) = {31, 29, 30};
Circle(22) = {30, 29, 31};
// 10:th circle 
Circle(23) = {34, 32, 33};
Circle(24) = {33, 32, 34};

// Loop over the lines to make them into one line
Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 6};
Line Loop(3) = {7, 8};
Line Loop(4) = {9, 10};
Line Loop(5) = {11, 12};
Line Loop(6) = {13, 14};
Line Loop(7) = {15, 16};
Line Loop(8) = {17, 18};
Line Loop(9) = {19, 20};
Line Loop(10) = {21, 22};
Line Loop(11) = {23, 24};

// Make the different surfaces 
Plane Surface(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};
Plane Surface(9) = {9};
Plane Surface(10) = {10};
Plane Surface(11) = {11};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Line(5) = {5, 6};
Physical Line(6) = {7, 8};
Physical Line(7) = {9, 10};
Physical Line(8) = {11, 12};
Physical Line(9) = {13, 14};
Physical Line(10) = {15, 16};
Physical Line(11) = {17, 18};
Physical Line(12) = {19, 20};
Physical Line(13) = {21, 22};
Physical Line(14) = {23, 24};

// Make the surface to a physical surface
Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};
Physical Surface(7) = {7};
Physical Surface(8) = {8};
Physical Surface(9) = {9};
Physical Surface(10) = {10};
Physical Surface(11) = {11};
