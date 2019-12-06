lc = 0.1;
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
	
	// Diff
	phiDiff = 0.20;
	
	//Add bottom point for the bottom part of the circle 
	phiVal = Atan2(y3, x3);
	theta = Pi / 2;
	x9 = R*Sin(theta)*Cos(phiVal); y9 = R*Sin(theta)*Sin(phiVal); z9 = R*Cos(theta);
	Point(pBot3) = {x9, y9, z9, lc}; 
	phiVal = Atan2(y1, x1);
	x10 = R*Sin(thetaHalf)*Cos(phiVal); y10 = R*Sin(thetaHalf)*Sin(phiVal); z10 = R*Cos(thetaHalf);

	// Add points for making trail to the bottom and upper segment for the right side point 
	phiVal = Atan2(y4, x4);
	theta = Pi / 2;
	x10 = R*Sin(theta)*Cos(phiVal); y10 = R*Sin(theta)*Sin(phiVal); z10 = R*Cos(theta);
	Point(pBot4) = {x10, y10, z10, lc};

	// For the right side point add small delta to make the mesh "whole"	
	phiVal = Atan2(y4, x4) + phiDiff;
	x10 = R*Sin(theta)*Cos(phiVal); y10 = R*Sin(theta)*Sin(phiVal); z10 = R*Cos(theta);
	Point(pBot5) = {x10, y10, z10, lc};
	// For this point add side point to build the trail
	theta = Acos(z4 / R); 
	x10 = R*Sin(theta)*Cos(phiVal); y10 = R*Sin(theta)*Sin(phiVal); z10 = R*Cos(theta);
	Point(pSide3) = {x10, y10, z10, lc};
	Point(pSide4) = {0, 0, z10, lc};
	theta = Acos(z1 / R);
	x10 = R*Sin(theta)*Cos(phiVal); y10 = R*Sin(theta)*Sin(phiVal); z10 = R*Cos(theta);
	Point(pSideUp2) = {x10, y10, z10, lc};
	// Add the extra side line 
	phiVal = Atan2(y4, x4) + 2*phiDiff;
	theta = Pi / 2;
	x10 = R*Sin(theta)*Cos(phiVal); y10 = R*Sin(theta)*Sin(phiVal); z10 = R*Cos(theta);
	Point(pSideRight1) = {x10, y10, z10, lc};
	theta = Acos(z4 / R);
	x10 = R*Sin(theta)*Cos(phiVal); y10 = R*Sin(theta)*Sin(phiVal); z10 = R*Cos(theta);
	Point(pSideRight2) = {x10, y10, z10, lc};
	 

	// Add point for making trail to the bottom for the left side point 
	phiVal = Atan2(y2, x2);
	theta = Pi / 2;
	x12 = R*Sin(theta)*Cos(phiVal); y12 = R*Sin(theta)*Sin(phiVal); z12 = R*Cos(theta);
	Point(pBot2) = {x12, y12, z12, lc}; 

	// For the left side add a small hole to move the mesh 
	phiVal = Atan2(y2, x2) - phiDiff;
	x12 = R*Sin(theta)*Cos(phiVal); y12 = R*Sin(theta)*Sin(phiVal); z12 = R*Cos(theta);
	Point(pBot1) = {x12, y12, z12, lc}; 
	theta = Acos(z2 / R); 
	x10 = R*Sin(theta)*Cos(phiVal); y10 = R*Sin(theta)*Sin(phiVal); z10 = R*Cos(theta);
	Point(pSide1) = {x10, y10, z10, lc};
	Point(pSide2) = {0, 0, z10, lc};
	theta = Acos(z1 / R);
	x10 = R*Sin(theta)*Cos(phiVal); y10 = R*Sin(theta)*Sin(phiVal); z10 = R*Cos(theta);
	Point(pSideUp1) = {x10, y10, z10, lc};
	// Point for making circle arcs for the upper side points 
	Point(pSideUp3) = {0, 0, z10, lc};
	// Point for the left support line 
	phiVal = Atan2(y2, x2) - 2* phiDiff;
	theta = Pi / 2;
	x10 = R*Sin(theta)*Cos(phiVal); y10 = R*Sin(theta)*Sin(phiVal); z10 = R*Cos(theta);
	Point(pSideLeft1) = {x10, y10, z10, lc};
	theta = Acos(z2 / R);
	x10 = R*Sin(theta)*Cos(phiVal); y10 = R*Sin(theta)*Sin(phiVal); z10 = R*Cos(theta);
	Point(pSideLeft2) = {x10, y10, z10, lc};

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

	// Add the path to the bottom for the circles 
	l9 = newl + 20;
	Circle(l9) = {p3, 1, pBot3};	
	l10 = newl + 20;
	Circle(l10) = {p4, 1, pBot4}; 
	l11 = newl + 20;
	Circle(l11) = {p2, 1, pBot2};

	// Circle arcs for the side points 
	l12 = newl;
	Circle(l12) = {pSide1, pSide2, p2};
	l13 = newl;
	Circle(l13) = {pSide3, pSide4, p4};
	l14 = newl;
	Circle(l14) = {pBot5, 1, pSide3};
	l15 = newl;
	Circle(l15) = {pBot1, 1, pSide1};
	l16 = newl;
	Circle(l16) = {pSide1, 1, pSideUp1};
	l17 = newl;
	Circle(l17) = {pSide3, 1, pSideUp2};
	l18 = newl;
	Circle(l18) = {pSideUp1, pSideUp3, p1};
	l19 = newl;
	Circle(l19) = {p1, pSideUp3, pSideUp2};

	// Circle arcs to the top-points 
	l20 = newl;
	Circle(l20) = {pSideUp1, 1, 6};
	l21 = newl;
	Circle(l21) = {p1, 1, 6};
	l22 = newl;
	Circle(l22) = {pSideUp2, 1, 6};

	// Circles for the extra side-points left
	l23 = newl;
	Circle(l23) = {pSideLeft1, 1, pSideLeft2};
	l24 = newl;
	Circle(l24) = {pSideLeft2, pSideUp3, pSide1};
	l25 = newl;
	Circle(l25) = {pSideLeft2, 1, 6};

	// Extra side point right 
	l26 = newl;
	Circle(l26) = {pSideRight1, 1, pSideRight2};
	l27 = newl;
	Circle(l27) = {pSideRight2, 1, 6};
	l28 = newl;
	Circle(l28) = {pSide3, pSide4, pSideRight2};

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

	// Fix arcs for the bottom parts 
	l29 = newl;
	Circle(l29) = {pSideLeft1, 1, pBot1};
	l30 = newl;
	Circle(l30) = {pBot1, 1, pBot2};
	l31 = newl;
	Circle(l31) = {pBot2, 1, pBot3};
	l32 = newl;
	Circle(l32) = {pBot3, 1, pBot4};
	l33 = newl;
	Circle(l33) = {pBot4, 1, pBot5};
	l34 = newl;
	Circle(l34) = {pBot5, 1, pSideRight1};

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
Point(7) = {R*Sin(thetaHalf)*Cos(0), R*Sin(thetaHalf)*Sin(0), R*Cos(thetaHalf), lc};
Point(8) = {R*Sin(thetaHalf)*Cos(Pi/2), R*Sin(thetaHalf)*Sin(Pi/2), R*Cos(thetaHalf), lc};
Point(9) = {R*Sin(thetaHalf)*Cos(Pi), R*Sin(thetaHalf)*Sin(Pi), R*Cos(thetaHalf), lc};
Point(10) = {R*Sin(thetaHalf)*Cos(3*Pi/2), R*Sin(thetaHalf)*Sin(3*Pi/2), R*Cos(thetaHalf), lc};

// Add a hole at the top
thetaHole = Acos(1 - 2*0.01);


// Add Hole 1
xRot = Pi / (2.9); yRot = Pi / 4.4; zRot = Pi / (2.50);
si = 101; psi = 2; 
pBot1 = 100; pBot2 = 101; pBot3 = 102; pBot4 = 103; pBot5 = 104;
pSide1 = 105; pSide2 = 106; pSide3 = 107; pSide4 = 108; pUp5 = 109;
pSideUp1 = 110; pSideUp2 = 111; pSideUp3 = 112;
pSideLeft1 = 113; pSideLeft2 = 114; pSideRight1 = 115; pSideRight2 = 116;
Call CreateHole;
pBot11 = pBot1; pBot12 = pBot2; pBot13 = pBot3; pBot14 = pBot4; pBot15 = pBot5;


// Building the bottom arc
Circle(1) = {3, 1, pSideLeft1};
Circle(2) = {pSideRight1, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Circle(5) = {2, 1, 3};

// Making arcs that go to the top 
Circle(30) = {2, 1, 6};
Circle(31) = {3, 1, 6};
Circle(32) = {4, 1, 6};
Circle(33) = {5, 1, 6};

// Fixing the curve-loops for the surfaces (via GUI)
Curve Loop(293) = {31, -285, -283, -1};
Curve Loop(294) = {283, 284, -275, -293};
Curve Loop(295) = {275, 272, 271, -294};
Curve Loop(296) = {82, 229, -295, -271};
Curve Loop(297) = {103, 250, -296, -229};
Curve Loop(298) = {250, 297, 274, 273};
Curve Loop(299) = {285, -280, -276, -284};
Curve Loop(300) = {278, 61, -272, 276};
Curve Loop(301) = {279, -277, 273, 124};
Curve Loop(302) = {278, 281, -280};
Curve Loop(303) = {279, 282, -281};
Curve Loop(304) = {287, -32, -2, 286};
Curve Loop(305) = {32, -33, -3};
Curve Loop(306) = {33, -30, -4};
Curve Loop(307) = {31, -30, 5};
Curve Loop(308) = {288, 287, -282, -277};
Curve Loop(309) = {298, 286, -288, -274};

// Fixing the different surfaces (via GUI)
Surface(293) = {293};
Surface(294) = {294};
Surface(295) = {295};
Surface(296) = {296};
Surface(297) = {297};
Surface(298) = {298};
Surface(299) = {299};
Surface(302) = {302};
Surface(303) = {303};
Surface(300) = {300};
Surface(301) = {301};
Surface(304) = {304};
Surface(305) = {305};
Surface(306) = {306};
Surface(307) = {307};
Surface(308) = {308};
Surface(309) = {309};

Physical Surface(1) = {293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309};


