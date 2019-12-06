lc = 0.01;
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
	x9 = R*Sin(theta)*Cos(phiVal); y9 = R*Sin(theta)*Sin(phiVal); z9 = R*Cos(theta);
	Point(pBot3) = {x9, y9, z9, lc}; 

	// Add points for making trail to the bottom and upper segment for the right side point 
	phiVal = Atan2(y4, x4);
	theta = Pi / 2;
	x10 = R*Sin(theta)*Cos(phiVal); y10 = R*Sin(theta)*Sin(phiVal); z10 = R*Cos(theta);
	Point(pBot4) = {x10, y10, z10, lc};

	// For the right side point add small delta to make the mesh "whole"	
	phiVal = Atan2(y4, x4) + 0.03;
	x10 = R*Sin(theta)*Cos(phiVal); y10 = R*Sin(theta)*Sin(phiVal); z10 = R*Cos(theta);
	Point(pBot5) = {x10, y10, z10, lc};

	// Add point for making trail to the bottom for the left side point 
	phiVal = Atan2(y2, x2);
	theta = Pi / 2;
	x12 = R*Sin(theta)*Cos(phiVal); y12 = R*Sin(theta)*Sin(phiVal); z12 = R*Cos(theta);
	Point(pBot2) = {x12, y12, z12, lc}; 

	// For the left side add a small hole to move the mesh 
	phiVal = Atan2(y2, x2) - 0.03;
	x12 = R*Sin(theta)*Cos(phiVal); y12 = R*Sin(theta)*Sin(phiVal); z12 = R*Cos(theta);
	Point(pBot1) = {x12, y12, z12, lc}; 

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

	// Add paths to the upper for the right side segment 
	l12 = newl + 20;
	Circle(l12) = {p1, 1, 6};
	l13 = newl + 20;
	Circle(l13) = {p4, 1, 6};
	l14 = newl + 20;
	Circle(l14) = {p2, 1, 6};

	// Add paths for the side regions 
	l15 = newl + 20;
	Circle(l15) = {pBot5, 1, 6};
	l16 = newl + 20;
	Circle(l16) = {pBot1, 1, 6};

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
thetaHole = Acos(1 - 2*0.01);

// Add Hole 1
xRot = Pi / (3.9); yRot = Pi / 2.60; zRot = Pi / (3*1.20);
si = 101; psi = 2; 
pBot1 = 100; pBot2 = 101; pBot3 = 102; pBot4 = 103; pBot5 = 104;
Call CreateHole;
pBot11 = pBot1; pBot12 = pBot2; pBot13 = pBot3; pBot14 = pBot4; pBot15 = pBot5;

// Add hole 2
xRot = Pi / (2.9); yRot = Pi / 4.4; zRot = Pi / (2.50);
si = 102; psi = 3; 
pBot1 = 200; pBot2 = 201; pBot3 = 202; pBot4 = 203; pBot5 = 204;
Call CreateHole;
pBot21 = pBot1; pBot22 = pBot2; pBot23 = pBot3; pBot24 = pBot4; pBot25 = pBot5;

// Add hole 3 
xRot = Pi / (2.45); yRot = Pi / 12.5; zRot = Pi / (2.25);
si = 103; psi = 4; 
pBot1 = 300; pBot2 = 301; pBot3 = 302; pBot4 = 303; pBot5 = 304;
Call CreateHole;
pBot31 = pBot1; pBot32 = pBot2; pBot33 = pBot3; pBot34 = pBot4; pBot35 = pBot5;

// Add hole 4 
xRot = Pi / (1.39); yRot = Pi / 1.65; zRot = Pi / (1.40);
si = 104; psi = 5; 
pBot1 = 400; pBot2 = 401; pBot3 = 402; pBot4 = 403; pBot5 = 404;
Call CreateHole;
pBot41 = pBot1; pBot42 = pBot2; pBot43 = pBot3; pBot44 = pBot4; pBot45 = pBot5;

// Add hole 5
xRot = Pi / (2.32); yRot = -Pi / 12.5; zRot = Pi / (2.0);
si = 105; psi = 6; 
pBot1 = 500; pBot2 = 501; pBot3 = 502; pBot4 = 503; pBot5 = 504;
Call CreateHole;
pBot51 = pBot1; pBot52 = pBot2; pBot53 = pBot3; pBot54 = pBot4; pBot55 = pBot5;


// Add the bottom circle arc 
Circle(1) = {5, 1, 2};
Circle(2) = {3, 1, pBot31};
Circle(3) = {pBot31, 1, pBot32};
Circle(4) = {pBot32, 1, pBot33};
Circle(5) = {pBot33, 1, pBot34};
Circle(6) = {pBot34, 1, pBot35};
Circle(7) = {pBot35, 1, pBot21};
Circle(8) = {pBot21, 1, pBot22};
Circle(9) = {pBot22, 1, pBot23};
Circle(10) = {pBot23, 1, pBot24};
Circle(11) = {pBot24, 1, pBot25};
Circle(12) = {pBot25, 1, pBot11};
Circle(13) = {pBot11, 1, pBot12};
Circle(14) = {pBot12, 1, pBot13};
Circle(15) = {pBot13, 1, pBot14};
Circle(16) = {pBot14, 1, pBot15};
Circle(17) = {pBot15, 1, 4};
Circle(18) = {4, 1, pBot41};
Circle(19) = {pBot41, 1, pBot42};
Circle(20) = {pBot42, 1, pBot43};
Circle(21) = {pBot43, 1, pBot44};
Circle(22) = {pBot44, 1, pBot45};
Circle(23) = {pBot45, 1, 5};
Circle(24) = {2, 1, pBot51};
Circle(25) = {pBot51, 1, pBot52};
Circle(26) = {pBot52, 1, pBot53};
Circle(27) = {pBot53, 1, pBot54};
Circle(28) = {pBot54, 1, pBot55};
Circle(29) = {pBot55, 1, 3};

// Making arcs that go to the top 
Circle(30) = {2, 1, 6};
Circle(31) = {3, 1, 6};
Circle(32) = {4, 1, 6};
Circle(33) = {5, 1, 6};

// Make arcs that go to the bottom 
p1 = newp;
Point(p1) = {0, 0, -R, lc};
Circle(34) = {pBot11, 1, p1};
Circle(35) = {pBot12, 1, p1};
Circle(36) = {pBot13, 1, p1};
Circle(37) = {pBot14, 1, p1};
Circle(38) = {pBot15, 1, p1};
Circle(39) = {pBot21, 1, p1};
Circle(40) = {pBot22, 1, p1};
Circle(41) = {pBot23, 1, p1};
Circle(42) = {pBot24, 1, p1};
Circle(43) = {pBot25, 1, p1};
Circle(44) = {pBot31, 1, p1};
Circle(45) = {pBot32, 1, p1};
Circle(46) = {pBot33, 1, p1};
Circle(47) = {pBot34, 1, p1};
Circle(48) = {pBot35, 1, p1};
Circle(49) = {pBot41, 1, p1};
Circle(50) = {pBot42, 1, p1};
Circle(51) = {pBot43, 1, p1};
Circle(52) = {pBot44, 1, p1};
Circle(53) = {pBot45, 1, p1};
Circle(54) = {pBot51, 1, p1};
Circle(55) = {pBot52, 1, p1};
Circle(56) = {pBot53, 1, p1};
Circle(57) = {pBot54, 1, p1};
Circle(58) = {2, 1, p1};
Circle(59) = {3, 1, p1};
Circle(60) = {4, 1, p1};
p2 = newp;
Circle(p2) = {5, 1, p1};

// The curve loops for the surface 
Curve Loop(1901) = {1896, -1854, 1791, -25};
Curve Loop(1902) = {26, -1749, -1602, 1791};
Curve Loop(1903) = {1581, 1854, -1812};
Curve Loop(1904) = {1644, 1812, -1833};
Curve Loop(1905) = {1623, 1770, -27, -1749};
Curve Loop(1906) = {1875, -1833, 1770, 28};
Curve Loop(1907) = {1875, -31, -29};
Curve Loop(1908) = {2, 1136, -31};
Curve Loop(1909) = {3, -1031, 1094, -1136};
Curve Loop(1910) = {4, -989, -842, 1031};
Curve Loop(1911) = {1094, -1052, 821};
Curve Loop(1912) = {884, 1052, -1073};
Curve Loop(1913) = {863, 1010, -5, -989};
Curve Loop(1914) = {1010, 6, 1115, -1073};
Curve Loop(1915) = {1115, -756, -7};
Curve Loop(1916) = {8, -651, 714, -756};
Curve Loop(1917) = {714, -672, 441};
Curve Loop(1918) = {651, 9, -609, -462};
Curve Loop(1919) = {483, 630, -10, -609};
Curve Loop(1920) = {672, -693, 504};
Curve Loop(1921) = {693, -735, -11, -630};
Curve Loop(1922) = {12, 376, -735};
Curve Loop(1923) = {13, -271, 334, -376};
Curve Loop(1924) = {82, 229, -14, -271};
Curve Loop(1925) = {229, 15, -250, -103};
Curve Loop(1926) = {124, 292, -313};
Curve Loop(1927) = {313, -355, -16, -250};
Curve Loop(1928) = {17, 32, -355};
Curve Loop(1929) = {1516, -32, 18};
Curve Loop(1930) = {19, -1411, 1474, -1516};
Curve Loop(1931) = {1201, 1474, -1432};
Curve Loop(1932) = {1222, 1369, -20, -1411};
Curve Loop(1933) = {1243, 1390, -21, -1369};
Curve Loop(1934) = {1432, -1453, 1264};
Curve Loop(1935) = {22, 1495, -1453, 1390};
Curve Loop(1936) = {23, 33, -1495};
Curve Loop(1937) = {30, -33, 1};
Curve Loop(1938) = {1896, -30, 24};
Curve Loop(1939) = {334, -292, 61};
Curve Loop(1940) = {53, -52, 22};
Curve Loop(1941) = {21, 52, -51};
Curve Loop(1942) = {20, 51, -50};
Curve Loop(1943) = {50, -49, 19};
Curve Loop(1944) = {49, -60, 18};
Curve Loop(1945) = {60, -38, 17};
Curve Loop(1946) = {38, -37, 16};
Curve Loop(1947) = {37, -36, 15};
Curve Loop(1948) = {36, -35, 14};
Curve Loop(1949) = {35, -34, 13};
Curve Loop(1950) = {34, -43, 12};
Curve Loop(1951) = {43, -42, 11};
Curve Loop(1952) = {42, -41, 10};
Curve Loop(1953) = {41, -40, 9};
Curve Loop(1954) = {40, -39, 8};
Curve Loop(1955) = {39, -48, 7};
Curve Loop(1956) = {48, -47, 6};
Curve Loop(1957) = {47, -46, 5};
Curve Loop(1958) = {46, -45, 4};
Curve Loop(1959) = {45, -44, 3};
Curve Loop(1960) = {44, -59, 2};
Curve Loop(1961) = {59, -57, 28, 29};
Curve Loop(1962) = {57, -56, 27};
Curve Loop(1963) = {56, -55, 26};
Curve Loop(1964) = {55, -54, 25};
Curve Loop(1965) = {24, 54, -58};
Curve Loop(1966) = {506, -58, -1};
Curve Loop(1967) = {506, -53, 23};

// The surfaces
Surface(1901) = {1901};
Surface(1902) = {1902};
Surface(1903) = {1903};
Surface(1904) = {1904};
Surface(1905) = {1905};
Surface(1906) = {1906};
Surface(1907) = {1907};
Surface(1908) = {1908};
Surface(1909) = {1909};
Surface(1910) = {1910};
Surface(1911) = {1911};
Surface(1912) = {1912};
Surface(1913) = {1913};
Surface(1914) = {1914};
Surface(1915) = {1915};
Surface(1916) = {1916};
Surface(1917) = {1917};
Surface(1918) = {1918};
Surface(1919) = {1919};
Surface(1920) = {1920};
Surface(1921) = {1921};
Surface(1922) = {1922};
Surface(1923) = {1923};
Surface(1924) = {1924};
Surface(1925) = {1925};
Surface(1926) = {1926};
Surface(1927) = {1927};
Surface(1928) = {1928};
Surface(1929) = {1929};
Surface(1930) = {1930};
Surface(1931) = {1931};
Surface(1932) = {1932};
Surface(1933) = {1933};
Surface(1934) = {1934};
Surface(1935) = {1935};
Surface(1936) = {1936};
Surface(1937) = {1937};
Surface(1938) = {1938};
Surface(1939) = {1939};
Surface(1940) = {1940};
Surface(1941) = {1941};
Surface(1942) = {1942};
Surface(1943) = {1943};
Surface(1944) = {1944};
Surface(1945) = {1945};
Surface(1946) = {1946};
Surface(1947) = {1947};
Surface(1948) = {1948};
Surface(1949) = {1949};
Surface(1950) = {1950};
Surface(1951) = {1951};
Surface(1952) = {1952};
Surface(1953) = {1953};
Surface(1954) = {1954};
Surface(1955) = {1955};
Surface(1956) = {1956};
Surface(1957) = {1957};
Surface(1958) = {1958};
Surface(1959) = {1959};
Surface(1960) = {1960};
Surface(1961) = {1961};
Surface(1962) = {1962};
Surface(1963) = {1963};
Surface(1964) = {1964};
Surface(1965) = {1965};
Surface(1966) = {1966};
Surface(1967) = {1967};

// Creating a physical surface 
Physical Surface(1) = {1901, 1902, 1903, 1904, 1905, 1906, 1907, 1908, 1909, 1910, 1911, 1912, 1913, 1914, 1915, 1916, 1917, 1918, 1919, 1920, 1921, 1922, 1923, 1924, 1925, 1926, 1927, 1928, 1929, 1930, 1931, 1932, 1933, 1934, 1935, 1936, 1937, 1938, 1939, 1940, 1941, 1942, 1943, 1944, 1945, 1946, 1947, 1948, 1949, 1950, 1951, 1952, 1953, 1954, 1955, 1956, 1957, 1958, 1959, 1960, 1961, 1962, 1963, 1964, 1965, 1966, 1967};
