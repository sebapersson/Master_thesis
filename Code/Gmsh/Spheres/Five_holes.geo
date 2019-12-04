lc = 0.05;
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
	x1 = R*Sin(thetaHole)*Cos(0); y1 = R*Sin(thetaHole)*Sin(0); z1 = R*Cos(thetaHole); 
	x2 = R*Sin(thetaHole)*Cos(Pi/2); y2 = R*Sin(thetaHole)*Sin(Pi/2); z2 = R*Cos(thetaHole); 
	x3 = R*Sin(thetaHole)*Cos(Pi); y3 = R*Sin(thetaHole)*Sin(Pi); z3 = R*Cos(thetaHole); 
	x4 = R*Sin(thetaHole)*Cos(3*Pi/2); y4 = R*Sin(thetaHole)*Sin(3*Pi/2); z4 = R*Cos(thetaHole); 
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


	//Add bottom point for the bottom part of the circle 
	phiVal = Atan2(y3, x3);
	theta = Pi / 2;
	x5 = R*Sin(theta)*Cos(phiVal); y5 = R*Sin(theta)*Sin(phiVal); z5 = R*Cos(theta);
	Point(pBot1) = {x5, y5, z5, lc}; 

	// Add points for making trail to the bottom for side point 
	phiVal = Atan2(y4, x4);
	theta = Pi / 2;
	x6 = R*Sin(theta)*Cos(phiVal); y6 = R*Sin(theta)*Sin(phiVal); z6 = R*Cos(theta);
	Point(pBot2) = {x6, y6, z6, lc}; 


	// Add circle on the tops, +10 to avoid conflicts with the bottom circle  
	l1 = newl + 20;
	Circle(l1) = {p1, pMean, p2};
	l2 = newl + 20;
	Circle(l2) = {p2, pMean, p3};
	l3 = newl + 20;
	Circle(l3) = {p3, pMean, p4};
	l4 = newl + 20;
	Circle(l4) = {p4, pMean, p1};
	// Add lines on the top 
	l5 = newl + 20;
	Circle(l5) = {p1, 1, pTop};
	l6 = newl + 20;
	Circle(l6) = {p2, 1, pTop};
	l7 = newl + 20;
	Circle(l7) = {p3, 1, pTop};
	l8 = newl + 20;
	Circle(l8) = {p4, 1, pTop};

	// Add the path to the bottom for the circles 
	l9 = newl + 20;
	Circle(l9) = {p3, 1, pBot1};	
	l10 = newl + 20;
	Circle(l10) = {p4, 1, pBot2}; 


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
Point(1) = {0, 0, 0, lc}; 
Point(2) = {R*Sin(theta)*Cos(0), R*Sin(theta)*Sin(0), 0, lc};
Point(3) = {R*Sin(theta)*Cos(Pi/2), R*Sin(theta)*Sin(Pi/2), 0, lc};
Point(4) = {R*Sin(theta)*Cos(Pi), R*Sin(theta)*Sin(Pi), 0, lc};
Point(5) = {R*Sin(theta)*Cos(3*Pi/2), R*Sin(theta)*Sin(3*Pi/2), 0, lc};
Point(6) = {0, 0, R, lc};

// Add a hole at the top
thetaHole = Pi / 9;
xRot = Pi / 4; yRot = Pi / 4; zRot = Pi / (3*1.35);

si = 101; psi = 2; pBot1 = 100; pBot2 = 1100;
Call CreateHole;

// Add the bottom circle arc 
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, pBot1};
Circle(3) = {pBot1, 1, pBot2};
Circle(4) = {pBot2, 1, 4};
Circle(5) = {4, 1, 5};
Circle(6) = {5, 1, 2};

// Making arcs that go to the top 
Circle(7) = {2, 1, 6};
Circle(8) = {3, 1, 6};
Circle(9) = {4, 1, 6};
Circle(10) = {5, 1, 6};

// Adding the second physical line 
Physical Line(1) = {1, 2, 3, 4, 5, 6};


