lc = 0.08;

// Outer rectangle points 
limitRec = 2.2;
Point(1) = {-limitRec, -limitRec, 0, lc};
Point(2) = {limitRec, -limitRec, 0, lc};
Point(3) = {limitRec, limitRec, 0, lc};
Point(4) = {-limitRec, limitRec, 0, lc};

// Lines between the points in the rectangle
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(1) = {1, 2, 3, 4};

Physical Line(1) = {1, 2, 3, 4};

// Creating a new circle macro 
Macro CreateHole
      // Three points define the circle
      p1 = newp; 
      Point(p1) = {x, y, 0, lc};
      p2 = newp; 
      Point(p2) = {x - r, y, 0, lc};
      p3 = newp;
      Point(p3) = {x + r, y, 0, lc};

      // Defining the circle curves
      c1 = newreg;
      Circle(c1) = {p3, p1, p2};
      c2 = newreg;
      Circle(c2) = {p2, p1, p3};

      // Making the curve lab 
      l1 = newreg;
      Line Loop(li) = {c1, c2};

      // Making the surface
      Plane Surface(s1) = {li};

      // Adding the physical line 
      Physical Line(li) = {c1, c2};
      Physical Surface(s1) = {s1};
      
Return
 
r = 0.25;

// Circle 1 
x = 0; y = 0; s1 = 2; li = 2;
Call CreateHole;

// Circle 2 
x = 0.4; y = 0.4; s1 = 3; li = 3;
Call CreateHole;

// Circle 3 
x = 0.55; y = -0.15; s1 = 4; li = 4;
Call CreateHole;

// Circle 4 
x = -0.15; y = -0.50; s1 = 5; li = 5;
Call CreateHole;

// Circle 5 
x = -0.55; y = 0.1; s1 = 6; li = 6;
Call CreateHole;

// Circle 6 
x = -0.15; y = 0.50; s1 = 7; li = 7;
Call CreateHole;

// Circle 7 
x = 0.40; y = -0.65; s1 = 8; li = 8;
Call CreateHole;

// Circle 8 
x = -0.70; y = -0.40; s1 = 9; li = 9;
Call CreateHole;

// Circle 9 
x = -0.50; y = -0.90; s1 = 10; li = 10;
Call CreateHole;

// Circle 10
x = 0.02; y = -1.0; s1 = 11; li = 11;
Call CreateHole;

// Circle 11
x = 0.20; y = 0.88; s1 = 12; li = 12;
Call CreateHole;

// Circle 12
x = -0.7; y = 0.60; s1 = 13; li = 13;
Call CreateHole;

// Circle 13
x = -1.1; y = 0.20; s1 = 14; li = 14;
Call CreateHole;

// Circle 14
x = 0.75; y = 0.85; s1 = 15; li = 15;
Call CreateHole;

// Circle 15
x = 0.95; y = 0.35; s1 = 16; li = 16;
Call CreateHole;

// Circle 16
x = -1.25; y = -0.30; s1 = 17; li = 17;
Call CreateHole;

// Circle 17 
x = -1.05; y = -0.85; s1 = 18; li = 18;
Call CreateHole;

// Circle 18 
x = 0.95; y = -0.55; s1 = 19; li = 19;
Call CreateHole;

// Circle 19 
x = -0.35; y = 1.0; s1 = 20; li = 20;
Call CreateHole;

// Circle 20 
x = -1.05; y = 1.0; s1 = 21; li = 21;
Call CreateHole;

Plane Surface(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
Physical Surface(1) = {1};

