// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/*************************************/
/* This example solves groundwater flow with a single well
 */
/*************************************/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <cmath>
#include <string>  
#include <iostream> 
#include <sstream>
// C++ includes
#include<math.h>
#include<iostream>
// dune-common includes
#include<dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include<dune/common/exceptions.hh> // We use exceptions
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/timer.hh>
// dune-geometry includes
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
// dune-grid includes
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/utility/structuredgridfactory.hh>
#include<dune/grid/io/file/vtk.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
// dune-pdelab includes
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/constraints/p0.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/function/callableadapter.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/variablefactories.hh>
#include<dune/pdelab/localoperator/l2.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestep.hh>
#include<dune/pdelab/newton/newton.hh>

// include stuff from dune-hydro
#include<dune/hydro/nonlineardiffusionfv.hh>
#include<dune/hydro/nonlineardiffusionfv_velocity.hh>

//GDAL includes. We require that GDAL is found
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()

// define the problem
template<typename Number>
class Model
{
  Dune::ParameterTree ptree;
  Number well;
  Number permeability;

public:
  //! Constructor without arg sets nonlinear term to zero
  Model (Dune::ParameterTree ptree_)
    : ptree(ptree_)
  {
    well = ptree.get("problem.well",(Number)0.0);
    permeability = ptree.get("problem.permeability",(Number)1.0);
  }

  //! storage term
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number storage (const E& e, const X& x) const
  {
    return 1.0;
  }

  //! bottom position in global coordinates (function b above)
  /**
     \param e      : codim 0 entity
     \param xlocal : local coordinate
   */
  template<typename E, typename X>
  Number bathymmetry (const E& e, const X& xlocal) const
  {
    return 0.0;
  }

  //! nonlinearity to be evaluated at a point on a face seperating two elements
  /**
     \param b_inside     : bottom position on inside element
     \param b_outside    : bottom position on outside element
     \param u_inside     : surface position on inside element
     \param u_outside    : surface position on outside element
     \param vn           : velocity in normal direction on face 
   */
  template<typename U, typename VN>
  Number phi (U b_inside, U b_outside, U u_inside, U u_outside, VN vn) const
  {
    // upwind evaluation of height
    Number u_upwind = (vn>=0) ? u_inside : u_outside;
    
    // evaluation of (u-b)
    Number bmax=std::max(b_inside,b_outside);
    Number h_upwind = std::max(u_upwind-bmax,0.0);
    
    return h_upwind;
  }

  template<typename U, typename VN>
  Number dphi_du_inside (U b_inside, U b_outside, U u_inside, U u_outside, VN vn) const
  {
    // check if we depend on u_inside
    if (vn<0) return 0.0;

    // upwind evaluation of height
    Number u_upwind = u_inside;

    // evaluation of (u-b)
    if (u_upwind-std::max(b_inside,b_outside)<0.0) return 0.0;
    
    return 1.0;
  }

  template<typename U, typename VN>
  Number dphi_du_outside (U b_inside, U b_outside, U u_inside, U u_outside, VN vn) const
  {
    // check if we depend on u_inside
    if (vn>=0) return 0.0;

    // upwind evaluation of height
    Number u_upwind = u_outside;

    // evaluation of (u-b)
    if (u_upwind-std::max(b_inside,b_outside)<0.0) return 0.0;

    return 1.0;
  }

  template<typename U, typename VN>
  Number dphi_dvn (U b_inside, U b_outside, U u_inside, U u_outside, VN vn) const
  {
    return 0.0;
  }

  //! permeability value at position in a cell
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number k (const E& e, const X& x) const
  {
    return permeability;
  }
  
  //! source/sink term
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number f (const E& e, const X& x) const
  {
    auto center = e.geometry().center();
    X domaincenter(50.0);
    domaincenter -= center;
    if (domaincenter.two_norm()<1e-5)
      return well;
    else
      return 0.0;
  }
  
  //! boundary condition type function (true = Dirichlet)
  /**
     \param i intersection (with boundary)
     \param x local coordinate on intersection
   */
  template<typename I, typename X>
  bool b (const I& i, const X& x) const
  {
    return false;
  }

  //! Dirichlet extension
  /** Provides Dirichlet boundary values as well as initial value!
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number g (const E& e, const X& x) const
  {
    return 10.0;
  }

  //! Neumann boundary condition
  /**
     \param i intersection (with boundary)
     \param x local coordinate on intersection
   */
  template<typename I, typename X>
  Number j (const I& i, const X& x) const
  {
    return 0.0;
  }

  //! set time for subsequent evaluation
  /**
     \param t value of time for subsequent evaluations
   */
  void setTime (double t)
  {}
};

//********************************************************************************
// driver function solving the problem
//********************************************************************************

template<typename GV, typename FEM>
void driverFV (const GV& gv, const FEM& fem, Dune::ParameterTree& ptree)
{
  // dimension and important types
  const int dim = GV::dimension;
  typedef double RF;                   // type for computations

  // Make grid function space
  typedef Dune::PDELab::P0ParallelConstraints CON;
  CON con;
  //typedef Dune::PDELab::NoConstraints CON;
  typedef Dune::PDELab::ISTL::VectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem,con);
  gfs.name("Vh");

  // A coefficient vector
  using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
  Z z(gfs); // initial value

  // Make a grid function out of it
  typedef Dune::PDELab::DiscreteGridFunction<GFS,Z> ZDGF;
  ZDGF zdgf(gfs,z);

  // make user functions and set initial time
  RF time = 0.0;
  using MODEL = Model<RF>;
  MODEL model(ptree);
  model.setTime(time);
  auto glambda = [&](const auto& e, const auto& x)
    {return model.g(e,x);};
  auto g = Dune::PDELab::
    makeInstationaryGridFunctionFromCallable(gv,glambda,model);
  auto blambda = [&](const auto& i, const auto& x)
    {return model.b(i,x);};
  auto b = Dune::PDELab::
    makeBoundaryConditionFromCallable(gv,blambda);
  auto bathymmetrylambda = [&](const auto& e, const auto& x)
    {return model.bathymmetry(e,x);};
  auto bathymmetrygf = Dune::PDELab::
    makeInstationaryGridFunctionFromCallable(gv,bathymmetrylambda,model);
  Z bathymmetry(gfs);
  Dune::PDELab::interpolate(bathymmetrygf,gfs,bathymmetry);
  ZDGF bdgf(gfs,bathymmetry);

  // Assemble constraints
  typedef typename GFS::template
    ConstraintsContainer<RF>::Type CC;
  CC cc; cc.clear();
  Dune::PDELab::constraints(b,gfs,cc); // assemble constraints
  std::cout << "constrained dofs=" << cc.size() << " of "
            << gfs.globalSize() << std::endl;

  // initialize simulation time,  the coefficient vector
  Dune::PDELab::interpolate(g,gfs,z);

  // Make instationary grid operator
  typedef NonlinearDiffusionFV<MODEL,GFS> LOP;
  LOP lop(model,gfs);
  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  MBE mbe(2*dim+1); // guess nonzeros per row
  //typedef Dune::PDELab::EmptyTransformation CC;
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,CC,CC> GO0;
  GO0 go0(gfs,cc,gfs,cc,lop,mbe);

  typedef FVL2<MODEL> TLOP;
  TLOP tlop(model);
  typedef Dune::PDELab::GridOperator<GFS,GFS,TLOP,MBE,RF,RF,RF,CC,CC> GO1;
  GO1 go1(gfs,cc,gfs,cc,tlop,mbe);

  typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;
  IGO igo(go0,go1);
  igo.divideMassTermByDeltaT();
  //igo.multiplySpatialTermByDeltaT();

  // Select a linear solver backend
  //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SSOR<IGO> LS;
  //typedef Dune::PDELab::ISTLBackend_SEQ_UMFPack LS;
  //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_ILU0 LS;
  //LS ls(1000,0);
  typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO> LS;
  LS ls (gfs,100,0);

  // solve nonlinear problem
  typedef Dune::PDELab::Newton<IGO,LS,Z> PDESOLVER;
  PDESOLVER pdesolver(igo,z,ls);
  pdesolver.setReassembleThreshold(0.0);
  pdesolver.setVerbosityLevel(2);
  pdesolver.setAbsoluteLimit(1e-10);
  pdesolver.setReduction(1e-8);
  pdesolver.setMinLinearReduction(1e-4);
  pdesolver.setMaxIterations(15);
  pdesolver.setLineSearchMaxIterations(10);
  // pdesolver.setLineSearchStrategy("hackbuschReuskenAcceptBest");
  pdesolver.setLineSearchStrategy("hackbuschReusken");

  // select and prepare time-stepping scheme
  Dune::PDELab::OneStepThetaParameter<RF> method1(1.0);
  Dune::PDELab::Alexander2Parameter<RF> method2;
  Dune::PDELab::Alexander3Parameter<RF> method3;
  int torder = ptree.get("method.torder",(int)1);
  Dune::PDELab::TimeSteppingParameterInterface<RF>*
    pmethod=&method1;
  if (torder==1) pmethod = &method1;
  if (torder==2) pmethod = &method2;
  if (torder==3) pmethod = &method3;
  if (torder<1||torder>3)
    std::cout<<"torder not in [1,3]"<<std::endl;
  typedef Dune::PDELab::OneStepMethod<RF,IGO,PDESOLVER,Z,Z> OSM;
  OSM  osm(*pmethod,igo,pdesolver);
  osm.setVerbosityLevel(2);

  // velocity postprocessing
  Z znew(z);
  ZDGF znewdgf(gfs,znew);
  typedef NonlinearDiffusionFVVelocity<MODEL,GFS,Z> VeloDGF;
  VeloDGF velodgf(model,gfs,znew,bathymmetry);

  // prepare VTK writer and write first file
  std::string filename=ptree.get("output.filename","output");
  struct stat st;
  if( stat( filename.c_str(), &st ) != 0 )
    {
      int stat = 0;
      stat = mkdir(filename.c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
      if( stat != 0 && stat != -1)
        std::cout << "Error: Cannot create directory "
                  << filename << std::endl;
    }
  typedef Dune::VTKWriter<GV> VTKWRITER;
  VTKWRITER vtkwriter(gv,Dune::VTK::conforming);
  typedef Dune::VTKSequenceWriter<GV> VTKSEQUENCEWRITER;
  VTKSEQUENCEWRITER vtkSequenceWriter(std::make_shared<VTKWRITER>(vtkwriter),filename,filename,"");
  typedef Dune::PDELab::VTKGridFunctionAdapter<ZDGF> VTKF;
  vtkSequenceWriter.addCellData(std::shared_ptr<VTKF>(new VTKF(znewdgf,"solution")));
  vtkSequenceWriter.addCellData(std::shared_ptr<VTKF>(new VTKF(bdgf,"bathymmetry")));
  typedef Dune::PDELab::VTKGridFunctionAdapter<VeloDGF> VeloVTKF;
  vtkSequenceWriter.addCellData(std::shared_ptr<VeloVTKF>(new VeloVTKF(velodgf,"velocity")));
  vtkSequenceWriter.write(time,Dune::VTK::appendedraw);

  // time loop
  RF T = ptree.get("problem.T",(RF)1.0);
  RF dt = ptree.get("method.dt",(RF)0.1);
  RF timestepmax =  ptree.get("method.dtmax",dt);
  int every = ptree.get("output.every",(int)1);
  int step=0;
  int increased_step = -1000;
  int decreased_step = -1000;
  double eps=1e-4;
  auto factor = 1.0/std::sqrt(2.0);
  while (time<T-1e-8)
    {
      // assemble constraints for new time step
      model.setTime(time+dt);

      // do time step
      try {
        znew = z;
        osm.apply(time,dt,z,znew);
      }
      catch (Dune::PDELab::NewtonLineSearchError) {
        dt *= factor;
        if (dt<eps) throw;
        decreased_step = step;
        if (gv.comm().rank()==0)
          std::cout << "STEP CONTROL: " << " decreased in step " << step+1 << " newdt: " << dt << std::endl;
        continue;
      }
      catch (Dune::PDELab::NewtonNotConverged) {
        dt *= factor;
        if (dt<eps) throw;
        decreased_step = step;
        if (gv.comm().rank()==0)
          std::cout << "STEP CONTROL: " << " decreased in step " << step+1 << " newdt: " << dt << std::endl;
        continue;
      }
      catch (Dune::PDELab::NewtonLinearSolverError) {
        dt *= factor;
        if (dt<eps) throw;
        decreased_step = step;
        if (gv.comm().rank()==0)
          std::cout << "STEP CONTROL: " << " decreased in step " << step+1 << " newdt: " << dt << std::endl;
        continue;
      }
      catch (Dune::ISTLError) {
        dt *= factor;
        if (dt<eps) throw;
        decreased_step = step;
        if (gv.comm().rank()==0)
          std::cout << "STEP CONTROL: " << " decreased in step " << step+1 << " newdt: " << dt << std::endl;
        continue;
      }

      // analyze change in time step
      double max_change = 0.0;
      auto& indexset = gv.indexSet();
      for (const auto& e : elements(gv,Dune::Partitions::interior)) {
        auto i = indexset.index(e);
        auto bi = Dune::PDELab::Backend::native(bathymmetry)[i][0];
        auto ziold = Dune::PDELab::Backend::native(z)[i][0];
        auto zinew = Dune::PDELab::Backend::native(znew)[i][0];
        max_change = std::max(max_change,std::abs(ziold-zinew));
      }
      gv.comm().max(max_change);
      if (gv.comm().rank()==0)
        std::cout << "MAXCHANGE " << max_change << " " << time+dt << " " << dt << std::endl;
        
      // accept time step
      velodgf.update();
      time+=dt;
      step++;
      z = znew;

      // determine maximum velicity value in any edge
      auto maxv = maxvelocity(velodgf);
      if (gv.comm().rank()==0)
        std::cout << "MAXVELOCITY " << maxv << std::endl;
      
      // output to VTK file
      if (step%every==0)
        {
          vtkSequenceWriter.write(time,Dune::VTK::appendedraw);
          if (gv.comm().rank()==0)
            std::cout << "WRITING VTK OUTPUT " << step << " " << time << " " << dt << std::endl;
        }
      
      // increase timestep
      if (step-decreased_step>3 && step-increased_step>3)
        {
          double newdt = std::min(timestepmax,std::sqrt(2.0)*dt);
          if (newdt>dt)
            {
              increased_step = step;
              if (gv.comm().rank()==0)
                std::cout << "STEP CONTROL: " << " increased in step " << step << " newdt: " << newdt << std::endl;
              dt = newdt;
            }
        }
    }
}


// main program: instantiates grid and calls driver
int main(int argc, char** argv)
{
  try{
    // Maybe initialize MPI
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    std::cout << "Hello World! This is dune-hydro." << std::endl;
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
      std::cout<<"I am rank "<<helper.rank()<<" of "<<helper.size()
               <<" processes!"<<std::endl;

    // register GDAL drivers once at the beginning
    GDALAllRegister();

    // open ini file
    Dune::ParameterTree ptree;
    Dune::ParameterTreeParser ptreeparser;
    ptreeparser.readINITree("firstexample.ini",ptree);
    ptreeparser.readOptions(argc,argv,ptree);

    // read ini file
    const int refinement = ptree.get<int>("grid.refinement");

    // read params
    const int dim=2;
    Dune::FieldVector<double,dim> L;
    L[0] = ptree.get("grid.structured.LX",(double)1.0);
    L[1] = ptree.get("grid.structured.LY",(double)1.0);
    std::array<int,dim> N;
    N[0] = ptree.get("grid.structured.NX",(unsigned int) 10);
    N[1] = ptree.get("grid.structured.NY",(unsigned int) 10);

    // YaspGrid section
    typedef Dune::YaspGrid<dim> Grid;
    typedef Grid::ctype DF;
    std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(new Grid(L,N));

    // UG structured triangular mesh
    // using Grid = Dune::UGGrid<dim>;
    // std::shared_ptr<Grid> gridp;
    // Dune::StructuredGridFactory<Grid> factory;
    // Dune::FieldVector<double,dim> origin;
    // origin[0] = origin[1] = 0.0;
    // gridp = factory.createSimplexGrid(origin,L,N);

    gridp->globalRefine(refinement);
    typedef Grid::LeafGridView GV;
    GV gv=gridp->leafGridView();

    typedef Dune::PDELab::P0LocalFiniteElementMap<double,double,dim> FEM;
    FEM fem(Dune::GeometryTypes::cube(dim));
    driverFV(gv,fem,ptree);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
