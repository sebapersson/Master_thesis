lc = 0.2;
R = 2;

// Rotate x, y, z for a given set of coordinates 
Macro RotateX
	x_new = x;
	y_new = y * Cos(xRot) + z * Sin(xRot);
	z_new = -y*Sin(xRot) + z * Cos(xRot);
	x = x_new; y = y_new; z = z_new;
Return 

Macro RotateY
	x_new = x * Cos(yRot) - z * Sin(yRot);
	y_new = y;
	z_new = x * Sin(yRot) + z * Cos(yRot);
	x = x_new; y = y_new; z = z_new;
Return 

Macro RotateZ
	x_new = x * Cos(zRot) + y * Sin(zRot);
	y_new = -x * Sin(zRot) + y * Cos(zRot);
	z_new = z;
	x = x_new; y = y_new; z = z_new;
Return 

// Building the base of the sphere
theta = Pi / 2.0;
Point(1) = {0, 0, 0}; 
Point(2) = {R*Sin(theta)*Cos(0), R*Sin(theta)*Sin(0), 0};
Point(3) = {R*Sin(theta)*Cos(Pi/2), R*Sin(theta)*Sin(Pi/2), 0};
Point(4) = {R*Sin(theta)*Cos(Pi), R*Sin(theta)*Sin(Pi), 0};
Point(5) = {R*Sin(theta)*Cos(3*Pi/2), R*Sin(theta)*Sin(3*Pi/2), 0};

// Add the bottom circle arc 
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Add a hole at the top
theta = Pi / 9;
theta2 = Pi / 18;
xRot = Pi / 4; yRot = Pi / 4; zRot = Pi / 4;

// Adding the points 
x1 = R*Sin(theta)*Cos(0); y1 = R*Sin(theta)*Sin(0); z1 = R*Cos(theta); 
x2 = R*Sin(theta)*Cos(Pi/2); y2 = R*Sin(theta)*Sin(Pi/2); z2 = R*Cos(theta); 
x3 = R*Sin(theta)*Cos(Pi); y3 = R*Sin(theta)*Sin(Pi); z3 = R*Cos(theta); 
x4 = R*Sin(theta)*Cos(3*Pi/2); y4 = R*Sin(theta)*Sin(3*Pi/2); z4 = R*Cos(theta); 
// The next set of points for the B-spline 
x5 = R*Sin(theta2)*Cos(0); y5 = R*Sin(theta2)*Sin(0); z5 = R*Cos(theta2); 
x6 = R*Sin(theta2)*Cos(Pi/2); y6 = R*Sin(theta2)*Sin(Pi/2); z6 = R*Cos(theta2); 
x7 = R*Sin(theta2)*Cos(Pi); y7 = R*Sin(theta2)*Sin(Pi); z7 = R*Cos(theta2); 
x8 = R*Sin(theta2)*Cos(3*Pi/2); y8 = R*Sin(theta2)*Sin(3*Pi/2); z8 = R*Cos(theta2); 
xTop = 0; yTop = 0; zTop = R;


// Perform the rotations 
// Point 1
x = x1; y = y1; z = z1;
Call RotateX; Call RotateY; Call RotateZ;
x1 = x; y1 = y; z1 = z;

// Point 2
x = x2; y = y2; z = z2;
Call RotateX; Call RotateY; Call RotateZ;
x2 = x; y2 = y; z2 = z;

// Point 3
x = x3; y = y3; z = z3;
Call RotateX; Call RotateY; Call RotateZ;
x3 = x; y3 = y; z3 = z;

// Point 4
x = x4; y = y4; z = z4;
Call RotateX; Call RotateY; Call RotateZ;
x4 = x; y4 = y; z4 = z;

// Point 5
x = x5; y = y5; z = z5;
Call RotateX; Call RotateY; Call RotateZ;
x5 = x; y5 = y; z5 = z;

// Point 6
x = x6; y = y6; z = z6;
Call RotateX; Call RotateY; Call RotateZ;
x6 = x; y6 = y; z6 = z;

// Point 7
x = x7; y = y7; z = z7;
Call RotateX; Call RotateY; Call RotateZ;
x7 = x; y7 = y; z7 = z;

// Point 8
x = x8; y = y8; z = z8;
Call RotateX; Call RotateY; Call RotateZ;
x8 = x; y8 = y; z8 = z;

// The top-point 
x = xTop; y = yTop; z = zTop;
Call RotateX; Call RotateY; Call RotateZ;
xTop = x; yTop = y; zTop = z;

xMean = (x1 + x2 + x3 + x4) / 4; yMean = (y1 + y2 + y3 + y4) / 4; zMean = (z1 + z2 + z3 + z4) / 4;

// Add the points 
Point(6) = {xTop, yTop, zTop, lc};
Point(7) = {x1, y1, z1, lc};
Point(8) = {x2, y2, z2, lc};
Point(9) = {x3, y3, z3, lc};
Point(10) = {x4, y4, z4, lc};
Point(11) = {xMean, yMean, zMean, lc};
Point(12) = {x5, y5, z5, lc};
Point(13) = {x6, y6, z6, lc};
Point(14) = {x7, y7, z7, lc};
Point(15) = {x8, y8, z8, lc};


// Add circle on the tops 
Circle(5) = {7, 11, 8};
Circle(6) = {8, 11, 9};
Circle(7) = {9, 11, 10};
Circle(8) = {10, 11, 7};
// Add lines on the top 
BSpline(9) = {9, 14, 6};
BSpline(10) = {8, 13, 6};
BSpline(11) = {7, 12, 6};
BSpline(12) = {10, 15, 6};

// Add the surfaces 
Curve Loop(1) = {9, -12, -7};
Surface(1) = {1};
Curve Loop(2) = {11, -12, 8};
Surface(2) = {2};
Curve Loop(3) = {5, 10, -11};
Surface(3) = {3};
Curve Loop(4) = {9, -10, 6};
Surface(4) = {4};
Surface Loop(1) = {1, 2, 3, 4};

// The curve that goes around
Curve Loop(5) = {5, 6, 7, 8};

Compoind Surface(30) = {1, 2, 3, 4};
