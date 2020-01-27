lc = 0.04;

R = 2.5;
// Outer rectangle points 
Point(1) = {0, -R, 0, lc};
Point(2) = {R, 0, 0, lc};
Point(3) = {0, R, 0, lc};
Point(4) = {-R, 0 , 0, lc};
Point(1000) = {0, 0, 0, lc};

// Lines between the points in the rectangle
Circle(1) = {1, 1000, 2};
Circle(2) = {2, 1000, 3};
Circle(3) = {3, 1000, 4};
Circle(4) = {4, 1000, 1};
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
x = 0.70; y = -0.12; s1 = 2; li = 2;
Call CreateHole;

// Circle 2 
x = 0.4; y = 0.4; s1 = 3; li = 3;
Call CreateHole;

// Circle 3 
x = 0.38; y = -0.6; s1 = 4; li = 4;
Call CreateHole;

// Circle 4
x = -0.15; y = -0.65; s1 = 5; li = 5;
Call CreateHole;

// Circle 5
x = -0.40; y = -0.15; s1 = 6; li = 6;
Call CreateHole;

// Circle 6
x = -0.15; y = 0.35; s1 = 7; li = 7;
Call CreateHole;

// Circle 7
x = 1.05; y = 0.300; s1 = 8; li = 8;
Call CreateHole;

Plane Surface(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Physical Surface(1) = {1};



