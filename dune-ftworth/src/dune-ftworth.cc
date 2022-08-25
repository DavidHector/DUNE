// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// always include the config file
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
// C++ includes
#include<math.h>
#include<iostream>
#include<vector>
#include<set>
#include<sstream>
// dune-common includes
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/timer.hh>
#include <dune/common/exceptions.hh> // We use exceptions
// dune-istl includes
#include<dune/istl/preconditioners.hh>
#include<dune/istl/paamg/amg.hh>
#include<dune/istl/paamg/pinfo.hh>
#include <dune/istl/umfpack.hh>
#include <dune/istl/cholmod.hh>
// dune-geometry includes
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
// dune-grid includes
#include<dune/grid/onedgrid.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include<dune/alugrid/grid.hh>
#include<dune/alugrid/dgf.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
// dune-pdelab includes
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/function/callableadapter.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
// ftworth includes
#include"../dune/ftworth/partitioner.hh"
#include"../dune/ftworth/assemblewrapper.hh"
#include"../dune/ftworth/subdomainutilities.hh"
#include"../dune/ftworth/cdediscretizer.hh"
#include"../dune/ftworth/schwarz_preconditioners.hh"
#include"../dune/ftworth/multilevel_geneo_preconditioner.hh"
// includes from this directory
#include"islands.hh"

#define DIMENSION 2
#define DEGREE 1
#define USECG 1
//#define USEDG 1
//#define USEFV 1
//#define UGGRID 1
#define YASPGRID 1

int main(int argc, char** argv)
{
  try{
    // Maybe initialize MPI
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    std::cout << "Hello World! This is dune-ftworth." << std::endl;
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
      std::cout<<"I am rank "<<helper.rank()<<" of "<<helper.size()
               <<" processes!"<<std::endl;

    // Read parameters from ini file
    Dune::ParameterTree ptree;
    Dune::ParameterTreeParser ptreeparser;
    ptreeparser.readINITree("dune-ftworth.ini",ptree);
    ptreeparser.readOptions(argc,argv,ptree);
   
#if (HAVE_UG) // check if we have the UG grid manager

    //===================
    // grid generation and partitioning
    //===================
    
    // load a grid file and refine
#ifdef UGGRID
    std::string meshtype = ptree.get("grid.type","structured");
    using Grid = Dune::UGGrid<DIMENSION>;  
    std::shared_ptr<Grid> pgrid;
    if (meshtype=="structured")
      {
        Dune::FieldVector<double,DIMENSION> lower_left(0.0);
        Dune::FieldVector<double,DIMENSION> upper_right;
        upper_right[0] = ptree.get("grid.LX",(double)1.0);
        upper_right[1] = ptree.get("grid.LY",(double)1.0);
        if (DIMENSION>2)
          upper_right[2] = ptree.get("grid.LZ",(double)1.0);
        std::array<unsigned int,DIMENSION> cells;
        cells[0] = ptree.get("grid.NX",(int)10);
        cells[1] = ptree.get("grid.NY",(int)10);
        if (DIMENSION>2)
          cells[2] = ptree.get("grid.NZ",(int)10);
        pgrid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(lower_left,upper_right,cells);
      }
    else
      {
        std::string mshfile = ptree.get("grid.filename","circle.msh");
        pgrid = std::shared_ptr<Grid>(Dune::GmshReader<Grid>::read(mshfile,true,false));
      }
#endif
#endif
    
    // make a YaspGrid
#ifdef YASPGRID
    using Grid = Dune::YaspGrid<DIMENSION>;  
    Dune::FieldVector<double,DIMENSION> upper_right;
    upper_right[0] = ptree.get("grid.LX",(double)1.0);
    upper_right[1] = ptree.get("grid.LY",(double)1.0);
    if (DIMENSION>2)
      upper_right[2] = ptree.get("grid.LZ",(double)1.0);
    std::array<int,DIMENSION> cells;
    cells[0] = ptree.get("grid.NX",(int)10);
    cells[1] = ptree.get("grid.NY",(int)10);
    if (DIMENSION>2)
      cells[2] = ptree.get("grid.NZ",(int)10);
    std::shared_ptr<Grid> pgrid = std::shared_ptr<Grid>(new Grid(upper_right,cells));
#endif

    // refine grid
    int refine = ptree.get("grid.refine",(int)0);
    pgrid->globalRefine(ptree.get("grid.refine",refine));

    // extract the leaf grid view
    auto gv = pgrid->leafGridView();
    using GV = decltype(gv); // and get the type
    const int dim = GV::dimension;
    auto& indexset = gv.indexSet(); // to attach data to elements
    std::vector<size_t> subdomains = ptree.get<std::vector<size_t>>("geneo.subdomains");
    std::cout << indexset.size(0) << " elements" << std::endl;
    std::cout << indexset.size(dim) << " vertices" << std::endl;
    std::cout << subdomains[0] << " subdomains" << std::endl;

    //===================
    // set up a discretization scheme
    //===================

    // solve a Poisson problem
    using DF = Grid::ctype;
    using RF = double;

    // set up the solver class
    double scaling = ptree.get("islands.scaling",(double)1.0);
    using Problem = IslandsProblem<GV,RF>;
    Problem problem(scaling);
#ifdef USECG
#ifdef UGGRID
    using Discretization = CGConvectionDiffusionProblemSimplex<GV,Problem,DEGREE>;
#endif
#ifdef YASPGRID
    using Discretization = CGConvectionDiffusionProblemCube<GV,Problem,DEGREE>;
#endif
    bool hasSkeleton = false;
#endif
#ifdef USEDG
#ifdef UGGRID
    using Discretization = DGConvectionDiffusionProblemSimplex<GV,Problem,DEGREE>;
#endif
#ifdef YASPGRID
    using Discretization = DGConvectionDiffusionProblemCube<GV,Problem,DEGREE>;
#endif
    bool hasSkeleton = true;
#endif
#ifdef USEFV
#ifdef YASPGRID
    using Discretization = FVConvectionDiffusionProblemCube<GV,Problem>;
#endif
    bool hasSkeleton = true;
#endif
    Discretization disc(gv,problem,ptree);

    //===================
    // set up global matrix and solution vector
    //===================

    // set up linear system in PDELab components
    using VEC = typename Discretization::VEC;
    using ISTLV = Dune::PDELab::Backend::Native<VEC>;
    using MAT = typename Discretization::MAT;
    using ISTLM = Dune::PDELab::Backend::Native<MAT>;

    VEC x(*disc.gfs);
    ISTLV& istlx = Dune::PDELab::Backend::native(x);
    Dune::PDELab::interpolate(*disc.g,*disc.gfs,x);
    MAT A(*disc.go); // this makes the sparsity pattern
    A = 0.0;
    Dune::Timer timer;
    disc.go->jacobian(x,A);
    auto assemble_time = timer.elapsed();
    std::cout << "ASSEMBLE TIME=" << assemble_time << std::endl;
    ISTLM& istlA = Dune::PDELab::Backend::native(A);
    std::shared_ptr<ISTLM> pistlA = stackobject_to_shared_ptr(Dune::PDELab::Backend::native(istlA));

    // residual and correction
    VEC r(*disc.gfs);
    r = 0.0;
    disc.go->residual(x,r);
    ISTLV& istlr = Dune::PDELab::Backend::native(r);

    // eliminate colum entries to Dirichlet nodes for Cholmod
    VEC v(*disc.gfs);
    ISTLV& istlv = Dune::PDELab::Backend::native(v);
    v = 1.0;
    set_constrained_dofs(*disc.cc,0.0,v); // dirichlet dofs are zero, others are 1, so we can use it as a mask!
	for (size_t i=0; i<istlA.N(); i++)
	  if (istlv[i][0]==1.0)
	    {
	      // non dirchlet row; eliminate Dirichlet columns
	      auto cIt = istlA[i].begin();
	      auto cEndIt = istlA[i].end();
	      for (; cIt!=cEndIt; ++cIt)
            if (istlv[cIt.index()][0]==0.0)
              (*cIt) = 0.0;
	    }
    
    //===================
    // now set up solver and solve
    //==================
    
    using SOLVER = Dune::Cholmod<ISTLV>;
    //using UMFPACKSOLVER = Dune::UMFPack<ISTLM>;
    timer.reset();
    //SOLVER solver(istlA,false);
    SOLVER solver;
    solver.setMatrix(istlA);
    auto factorization_time = timer.elapsed();
    std::cout << "FACTORIZATION TIME=" << factorization_time << std::endl;
    Dune::InverseOperatorResult stat;
    for (size_t i=0; i<istlv.N(); i++)
      for (size_t j=0; j<istlv[i].N(); j++)
        istlv[i][j] = (rand() % 10000)/10000.0;
    set_constrained_dofs(*disc.cc,0.0,v); // dirichlet dofs are zero, others are 1, so we can use it as a mask!
    timer.reset();    
    solver.apply(istlv,istlr,stat);
    auto solver_time = timer.elapsed();
    std::cout << "SOLVER TIME=" << solver_time << std::endl;
    istlx -= istlv;

    std::cout << "COMBINED "
              << " dof=" << disc.gfs->globalSize()
              << " ass=" << assemble_time
              << " fac=" << factorization_time
              << " slv=" << solver_time
              << std::endl;
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
