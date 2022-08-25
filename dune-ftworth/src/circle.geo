// mesh width associated with points
Mesh.MshFileVersion = 2.2;
radius = 1.0;
meshsize = 0.03;

Point(1) = {0, 0, 0, meshsize};
Point(2) = {radius, 0, 0, meshsize};
Point(3) = {-radius, 0, 0, meshsize};

Circle(1) = {2,1,3};
Circle(2) = {3,1,2};

Curve Loop(100) = {1,2};  
Plane Surface(200) = {100};  
