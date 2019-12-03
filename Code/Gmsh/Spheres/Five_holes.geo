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

Macro CreateHole
	
	// Adding the points used to create a hole  
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

	// Calculate the mean for constructing the circle 
	xMean = (x1 + x2 + x3 + x4) / 4; yMean = (y1 + y2 + y3 + y4) / 4; zMean = (z1 + z2 + z3 + z4) / 4;

	// Add the points 
	p1 = newp;
	Point(p1) = {x1, y1, z1, lc};
	p2 = newp;
	Point(p2) = {x2, y2, z2, lc};
	p3 = newp;
	Point(p3) = {x3, y3, z3, lc};
	p4 = newp;
	Point(p4) = {x4, y4, z4, lc};
	pMean = newp;
	Point(pMean) = {xMean, yMean, zMean, lc};
	pTop = newp;
	Point(pTop) = {xTop, yTop, zTop, lc};
	p5 = newp;
	Point(p5) = {x5, y5, z5, lc};
	p6 = newp;
	Point(p6) = {x6, y6, z6, lc};
	p7 = newp;
	Point(p7) = {x7, y7, z7, lc};
	p8 = newp;
	Point(p8) = {x8, y8, z8, lc};

	// Add points for making a trail to the top 
	phiVal = Atan2(y1, x1);
	thetaVal = Acos(z1 / R);
	theta1 = thetaVal / 1.2;
	theta2 = thetaVal / 2;
	theta3 = thetaVal / 4;
	x9 = R*Sin(theta1)*Cos(phiVal); y9 = R*Sin(theta1)*Sin(phiVal); z9 = R*Cos(theta1); 
	x10 = R*Sin(theta2)*Cos(phiVal); y10 = R*Sin(theta2)*Sin(phiVal); z10 = R*Cos(theta2); 
	x11 = R*Sin(theta3)*Cos(phiVal); y11 = R*Sin(theta3)*Sin(phiVal); z11 = R*Cos(theta3);
	p9 = newp;
	Point(p9) = {x9, y9, z9}; 
	p10 = newp;
	Point(p10) = {x10, y10, z10}; 
	p11 = newp;
	Point(p11) = {x11, y11, z11}; 

	//Add points for making trail to the bottom 
	phiVal = Atan2(y3, x3);
	thetaVal = Acos(z3 / R);
	diff = Pi/2 - thetaVal;
	theta1 = (thetaVal + diff / 1.2);
	theta2 = (thetaVal + diff / 2);
	theta3 = (thetaVal + diff / 3);
	theta4 = Pi / 2;
	x12 = R*Sin(theta1)*Cos(phiVal); y12 = R*Sin(theta1)*Sin(phiVal); z12 = R*Cos(theta1); 
	x13 = R*Sin(theta2)*Cos(phiVal); y13 = R*Sin(theta2)*Sin(phiVal); z13 = R*Cos(theta2); 
	x14 = R*Sin(theta3)*Cos(phiVal); y14 = R*Sin(theta3)*Sin(phiVal); z14 = R*Cos(theta3);
	x15 = R*Sin(theta4)*Cos(phiVal); y15 = R*Sin(theta4)*Sin(phiVal); z15 = R*Cos(theta4);
	p12 = newp;
	Point(p12) = {x12, y12, z12}; 
	p13 = newp;
	Point(p13) = {x13, y13, z13}; 
	p14 = newp;
	Point(p14) = {x14, y14, z14};
	Point(pBot) = {x15, y15, z15}; 


	// Add circle on the tops, +1000 to avoid conflicts with the bottom circle  
	l1 = newl + 1000;
	Circle(l1) = {p1, pMean, p2};
	l2 = newl + 1000;
	Circle(l2) = {p2, pMean, p3};
	l3 = newl + 1000;
	Circle(l3) = {p3, pMean, p4};
	l4 = newl + 1000;
	Circle(l4) = {p4, pMean, p1};
	// Add lines on the top 
	l5 = newl + 1000;
	BSpline(l5) = {p1, p5, pTop};
	l6 = newl + 1000;
	BSpline(l6) = {p2, p6, pTop};
	l7 = newl + 1000;
	BSpline(l7) = {p3, p7, pTop};
	l8 = newl + 1000;
	BSpline(l8) = {p4, p8, pTop};

	// For the trail up to the top 
	l9 = newl + 1000; 
	BSpline(l9) = {p1, p9, p10, p11, 6};
	// For the trail to the bottom
	l10 = newl + 1000;
	BSpline(l10) = {p3, p12, p13, p14, pBot};

	// Add the surfaces 
	c1 = newreg; 
	Curve Loop(c1) = {l1, l6, -l5};
	Surface(c1) = {c1};
	c2 = newreg; 
	Curve Loop(c2) = {l2, l7, -l6};
	Surface(c2) = {c2};
	c3 = newreg; 
	Curve Loop(c3) = {l3, l8, -l7};
	Surface(c3) = {c3};
	c4 = newreg; 
	Curve Loop(c4) = {l4, l5, -l8};
	Surface(c4) = {c4};

	// For making the big surface in a later stage 
	Curve Loop(si) = {l1, l2, l3, l4};

	// Fix the physical entities
	Physical Surface(psi) = {c1, c2, c3, c4};	
	Physical Line(psi) = {l1, l2, l3, 4};

Return 

// Building the base of the sphere
theta = Pi / 2.0;
Point(1) = {0, 0, 0}; 
Point(2) = {R*Sin(theta)*Cos(0), R*Sin(theta)*Sin(0), 0};
Point(3) = {R*Sin(theta)*Cos(Pi/2), R*Sin(theta)*Sin(Pi/2), 0};
Point(4) = {R*Sin(theta)*Cos(Pi), R*Sin(theta)*Sin(Pi), 0};
Point(5) = {R*Sin(theta)*Cos(3*Pi/2), R*Sin(theta)*Sin(3*Pi/2), 0};
Point(6) = {0, 0, R};



// Add a hole at the top
theta = Pi / 9;
theta2 = Pi / 18;
xRot = Pi / 4; yRot = Pi / 4; zRot = Pi / 3;

si = 101; psi = 2; pBot = 100;
Call CreateHole;

// Add the bottom circle arc 
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, pBot};
Circle(3) = {pBot, 1, 4};
Circle(4) = {4, 1, 5};
Circle(5) = {5, 1, 2};

// Making arcs that go to the top 
Circle(6) = {2, 1, 6};
Circle(7) = {3, 1, 6};
Circle(8) = {4, 1, 6};
Circle(9) = {5, 1, 6};

// Adding the second physical line 
Physical Line(1) = {1, 2, 3, 4, 5};


