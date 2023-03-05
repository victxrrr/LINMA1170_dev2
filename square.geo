cl = 1.;

Point(1) = {0.,0.,0., cl};
Point(2) = {1.,0.,0., cl};
Point(3) = {1.,1.,0., cl};
Point(4) = {0.,1.,0., cl};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Surface(1) = {1};

Physical Curve("clamped", 1) = {4};
//+
Physical Surface("bulke", 1) = {1};
//+
//Physical Point("forcex", 7) = {2};
//+
Physical Line("forcey", 8) = {2};

