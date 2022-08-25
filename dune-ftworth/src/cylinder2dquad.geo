// Original file by Marian Piatkowski
// modified by Peter Bastian
Mesh.MshFileVersion = 2.2;

nPointsStandard = 17;
nPointsLong = 25;
nPointsBndLayer = 8;
nPointsRadial = 17;
progressionRadial = 0.85;
progressionBndLayer = 1.6;
progressionLong = 1.04;


Point(1) = {0     , 0.06, 0, 1.0};
Point(2) = {0.41, 0.06, 0, 1.0};
Point(3) = {1.41, 0.06, 0, 1.0};
Point(4) = {1.41, 0.35, 0, 1.0};
Point(5) = {0.41, 0.35, 0, 1.0};
Point(6) = {0     , 0.35, 0, 1.0};

Point(7) = {0.2, 0.2, 0, 1.0};
Point(8) = {0.2+0.05*Cos(Pi/4), 0.2+0.05*Sin(Pi/4), 0, 1.0};
Point(9) = {0.2+0.05*Cos(3*Pi/4), 0.2+0.05*Sin(3*Pi/4), 0, 1.0};
Point(10) = {0.2+0.05*Cos(5*Pi/4), 0.2+0.05*Sin(5*Pi/4), 0, 1.0};
Point(11) = {0.2+0.05*Cos(7*Pi/4), 0.2+0.05*Sin(7*Pi/4), 0, 1.0};

Point(12) = {0     , 0, 0, 1.0};
Point(13) = {0.41, 0, 0, 1.0};
Point(14) = {1.41, 0, 0, 1.0};
Point(15) = {1.41, 0.41, 0, 1.0};
Point(16) = {0.41, 0.41, 0, 1.0};
Point(17) = {0     , 0.41, 0, 1.0};

Line(1) = {1, 2}; // box out
Line(3) = {6, 5}; // box out
Line(5) = {1, 6}; // inlet
Line(6) = {2, 5}; // middle
Line(7) = {3, 4}; // outlet transversal

Line(2) = {2, 3}; // long edge box to end
Line(4) = {5, 4}; // long edge box to end

Circle(8) = {10, 7, 11}; // correp. to (1,2)=1
Circle(9) = {11, 7, 8};  // correp. to (2,5)=6
Circle(10) = {8, 7, 9};  // correp. to (5,6)=3
Circle(11) = {9, 7, 10}; // correp. to (6,1)=5

Line(12) = {1, 10};
Line(13) = {2, 11};
Line(14) = {5, 8};
Line(15) = {6, 9};

Line(16) = {12,1};
Line(17) = {12,13};
Line(18) = {13,2};
Line(19) = {13,14};
Line(20) = {14,3};

Line(21) = {15,4};
Line(22) = {16,15};
Line(23) = {16,5};
Line(24) = {16,17};
Line(25) = {17,6};


Transfinite Line {1} = nPointsStandard Using Progression 1;
Transfinite Line {3} = nPointsStandard Using Progression 1;
Transfinite Line {5} = nPointsStandard Using Progression 1;
Transfinite Line {6} = nPointsStandard Using Progression 1;
Transfinite Line {7} = nPointsStandard Using Progression 1;
Transfinite Line {8} = nPointsStandard Using Progression 1;
Transfinite Line {9} = nPointsStandard Using Progression 1;
Transfinite Line {10} = nPointsStandard Using Progression 1;
Transfinite Line {11} = nPointsStandard Using Progression 1;

Transfinite Line {2} = nPointsLong Using Progression progressionLong; // long part
Transfinite Line {4} = nPointsLong Using Progression progressionLong; // long part

Transfinite Line {12} = nPointsRadial Using Progression progressionRadial; //radial edges
Transfinite Line {13} = nPointsRadial Using Progression progressionRadial;
Transfinite Line {14} = nPointsRadial Using Progression progressionRadial;
Transfinite Line {15} = nPointsRadial Using Progression progressionRadial;

Transfinite Line {17} = nPointsStandard Using Progression 1;
Transfinite Line {24} = nPointsStandard Using Progression 1;

Transfinite Line {19} = nPointsLong Using Progression progressionLong;
Transfinite Line {22} = nPointsLong Using Progression progressionLong;

Transfinite Line {16} = nPointsBndLayer Using Progression progressionBndLayer;
Transfinite Line {18} = nPointsBndLayer Using Progression progressionBndLayer;
Transfinite Line {20} = nPointsBndLayer Using Progression progressionBndLayer;
Transfinite Line {21} = nPointsBndLayer Using Progression progressionBndLayer;
Transfinite Line {23} = nPointsBndLayer Using Progression progressionBndLayer;
Transfinite Line {25} = nPointsBndLayer Using Progression progressionBndLayer;

Line Loop(16) = {12, 8, -13, -1};
Plane Surface(17) = {16};
Line Loop(18) = {13, 9, -14, -6};
Plane Surface(19) = {18};
Line Loop(20) = {14, 10, -15, 3};
Plane Surface(21) = {20};
Line Loop(22) = {15, 11, -12, 5};
Plane Surface(23) = {22};
Line Loop(24) = {2, 7, -4, -6};
Plane Surface(25) = {24};

Line Loop(26) = {-16, 17, 18, -1};
Plane Surface(27) = {26};
Line Loop(28) = {-18, 19, 20, -2};
Plane Surface(29) = {28};
Line Loop(30) = {23, 4, -21, -22};
Plane Surface(31) = {30};
Line Loop(32) = {25, 3, -23, 24};
Plane Surface(33) = {32};

Transfinite Surface {17};
Transfinite Surface {19};
Transfinite Surface {21};
Transfinite Surface {23};
Transfinite Surface {25};
Transfinite Surface {27};
Transfinite Surface {29};
Transfinite Surface {31};
Transfinite Surface {33};


Recombine Surface {17};
Recombine Surface {19};
Recombine Surface {21};
Recombine Surface {23};
Recombine Surface {25};
Recombine Surface {27};
Recombine Surface {29};
Recombine Surface {31};
Recombine Surface {33};

//Extrude {0, 0, 0.41} {
//  Surface{23, 17, 19, 21, 25};
//  Layers{5};
//  Recombine;
//}
