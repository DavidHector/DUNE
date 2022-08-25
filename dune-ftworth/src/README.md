# Documentation for directory arc
This directory contains specific examples for testing the multilevel GENEO variants.

##The different examples

###testproblems.cc
Is a simple Poisson problem defined in Poisson.hh with homogeneous Dirichlet boundary conditions and right hand side 1 which will work on any geometry in any dimension. The file can be easily edited to solve different problems. Just
- change the dimension of the grid UGGrid<2> or UGGrid<3>
- use any of the classes
  * CGConvectionDiffusionProblemCube - conforming finite element method on cubes
  * CGConvectionDiffusionProblemSimplex - conforming finite element method on simplices
  * DGConvectionDiffusionProblemCube - DG finite elements on cubes
  * DGConvectionDiffusionProblemSimplex - DG finite elements on simplices
The classes take also the polynomial degree as a parameter.

###spe10.cc
Solves the SPE10 test problem on a fixed hexahedral mesh of 60x220x85 elements. It can be refined though.
Before you can use the SPE10 example you have to extract the archive file spe10data.tgz and move the file spe10_perm.dat
to the directory where the executable lies.

