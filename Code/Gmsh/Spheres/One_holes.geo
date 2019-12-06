lc = 0.1;
R = 2;


Macro CreateHole
	
	// Adding the points used to create a hole  
	x1 = R*Sin(thetaHole)*Cos(0); y1 = R*Sin(thetaHole)*Sin(0); z1 = R*Cos(thetaHole); 
	x2 = R*Sin(thetaHole)*Cos(Pi/2); y2 = R*Sin(thetaHole)*Sin(Pi/2); z2 = R*Cos(thetaHole); 
	x3 = R*Sin(thetaHole)*Cos(Pi); y3 = R*Sin(thetaHole)*Sin(Pi); z3 = R*Cos(thetaHole); 
	x4 = R*Sin(thetaHole)*Cos(3*Pi/2); y4 = R*Sin(thetaHole)*Sin(3*Pi/2); z4 = R*Cos(thetaHole); 
	xTop = 0; yTop = 0; zTop = R;

	// Calculate the mean for constructing the circle 
	xMean = (x1 + x2 + x3 + x4) / 4; yMean = (y1 + y2 + y3 + y4) / 4; zMean = (z1 + z2 + z3 + z4) / 4;

	// Add the points 
	p1 = newp;
	Point(p1) = {x1, y1, z1, lc};
	Printf("x1^2 + y1^2 + z1^2 = %.3f, p1 = %.3f\n", x1*x1+y1*y1+z1*z1, p1);
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
	

	// Add circle on the tops, +10 to avoid conflicts with the bottom circle  
	l1 = newl + 60;
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
thetaHalf = Pi / 3;
Point(1) = {0, 0, 0, lc}; 
Point(2) = {R*Sin(theta)*Cos(0), R*Sin(theta)*Sin(0), 0, lc};
Point(3) = {R*Sin(theta)*Cos(Pi/2), R*Sin(theta)*Sin(Pi/2), 0, lc};
Point(4) = {R*Sin(theta)*Cos(Pi), R*Sin(theta)*Sin(Pi), 0, lc};
Point(5) = {R*Sin(theta)*Cos(3*Pi/2), R*Sin(theta)*Sin(3*Pi/2), 0, lc};
Point(6) = {0, 0, R, lc};


// Add a hole at the top
thetaHole = Acos(1 - 2*0.01);


// Add Hole 1
xRot = Pi / (2.9); yRot = Pi / 4.4; zRot = Pi / (2.50);
si = 101; psi = 2; 
Call CreateHole;

p = newp;
Point(p) = {0, 0, -R, lc};





//+
Circle(209) = {5, 1, 2};
//+
Circle(210) = {2, 1, 3};
//+
Circle(211) = {3, 1, 4};
//+
Circle(212) = {4, 1, 5};
//+
Circle(213) = {4, 1, 9};
//+
Circle(214) = {5, 1, 10};
//+
Circle(215) = {2, 1, 7};
//+
Circle(216) = {3, 1, 8};
//+
Circle(217) = {3, 1, 13};
//+
Circle(218) = {2, 1, 13};
//+
Circle(219) = {5, 1, 13};
//+
Circle(220) = {4, 1, 13};
//+
Curve Loop(213) = {213, 103, -214, -212};
//+
Surface(213) = {213};
//+
Curve Loop(214) = {209, 215, -124, -214};
//+
Surface(214) = {214};
//+
Curve Loop(215) = {215, 61, -216, -210};
//+
Surface(215) = {215};
//+
Curve Loop(216) = {216, 82, -213, -211};
//+
Surface(216) = {216};
//+
Curve Loop(217) = {211, 220, -217};
//+
Surface(217) = {217};
//+
Curve Loop(218) = {219, -220, 212};
//+
Surface(218) = {218};
//+
Curve Loop(219) = {209, 218, -219};
//+
Surface(219) = {219};

//+
Curve Loop(220) = {210, 217, -218};
//+
Surface(220) = {220};

Physical Surface(1) = {213, 214, 215, 216, 217, 218, 219, 220};
