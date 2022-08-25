// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/*************************************/
/* This example solves surface flow on an inclined plane with transport
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
#include<dune/pdelab/finiteelementmap/qkdg.hh>
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
#include <dune/pdelab/localoperator/convectiondiffusiondg.hh>
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

//********************************************************************************
// Parameter class for flow part
//********************************************************************************

template<typename Number>
class Model
{
  Number eps,A,B;

public:
  //! Constructor without arg sets nonlinear term to zero
  Model () : eps(1e-3)
  {
    // regularization parameter
    A = 2.5/(2.0*std::sqrt(eps));
    B = 0.5/(2.0*pow(eps,2.5));
    std::cout << "REGULARIZATION gradient eps=" << eps << " A=" << A << " B=" << B << std::endl;
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
    auto x = e.geometry().center();
    return 5.0+0.01*x[0]+0.04*x[1]; // the inclined plane
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
    
    // regularization of gradient and evaluation
    Number gradterm =  (std::abs(vn)>=eps) ? 1.0/std::sqrt(std::abs(vn)) : (A-B*vn*vn);

    return h_upwind*gradterm;
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
    
    // regularization of gradient and evaluation
    Number gradterm =  (std::abs(vn)>=eps) ? 1.0/std::sqrt(std::abs(vn)) : (A-B*vn*vn);

    return gradterm;
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
    
    // regularization of gradient and evaluation
    Number gradterm =  (std::abs(vn)>=eps) ? 1.0/std::sqrt(std::abs(vn)) : (A-B*vn*vn);

    return gradterm;
  }

  template<typename U, typename VN>
  Number dphi_dvn (U b_inside, U b_outside, U u_inside, U u_outside, VN vn) const
  {
    // upwind evaluation of height
    Number u_upwind = (vn>=0) ? u_inside : u_outside;

    // evaluation of (u-b)
    Number bmax=std::max(b_inside,b_outside);
    Number h_upwind = std::max(u_upwind-bmax,0.0);
    if (h_upwind<0.0) h_upwind = 0.0;
    
    // regularization of gradient and evaluation
    Number sign = (vn>=0) ? 1.0 : -1.0;
    if (std::abs(vn)>=eps)
      return -h_upwind*0.5/(std::abs(vn)*std::sqrt(std::abs(vn)))*sign;
    else
      return -h_upwind*B*2.0*vn;
  }

  //! permeability value at position in a cell
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number k (const E& e, const X& x) const
  {
    return 1.0;
  }
  
  //! source/sink term
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number f (const E& e, const X& x) const
  {
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
    auto global = e.geometry().global(x);
    decltype(global) center;
    center[0] = 50.0;
    center[1] = 25.0;
    center -= global;
    Number rmax = 15.0;
    Number r = center.two_norm();
    return bathymmetry(e,x) + 2*(0.5+atan(2*(rmax-r))/M_PI);
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
// Parameter class for transport part
//********************************************************************************

template<typename GV, typename P0DGF, typename VDGF>
class ModelTransport
{
  using RF = typename VDGF::Traits::RangeFieldType;
  using BCType = typename Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type;
  using DiffusionTensor =  typename Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF>::PermTensorType;

  DiffusionTensor D;
  mutable RF time;
  mutable RF starttime;
  mutable RF dt;
  const P0DGF& bathymmetry;
  const P0DGF& elevationold;
  const P0DGF& elevationnew;
  const VDGF& velocity;
  RF residualheight;
  
public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  ModelTransport (const Dune::ParameterTree& ptree, const P0DGF& bathymmetry_,
                  const P0DGF& elevationold_, const P0DGF& elevationnew_, const VDGF& velocity_)
    : time(0.0), bathymmetry(bathymmetry_), elevationold(elevationold_), elevationnew(elevationnew_), velocity(velocity_)
  {
    auto myD = ptree.get<RF>("problem.D");
    residualheight = ptree.get<RF>("problem.residualheight");

    std::cout << "MODELT Constructor" << " D=" << D[0][0] << " residualheight=" << residualheight << std::endl;

    // diffusion tensor
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        D[i][j] = (i==j) ? myD : 0;
  }

  //! storage term
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  RF storage (const E& e, const X& x) const
  {
    auto cellcenterlocal = referenceElement(e.geometry()).position(0,0);

    // check time
    // if (!(std::abs(time-starttime)<1e-13) && !(std::abs(time-starttime-dt)<1e-13))
    //   {
    //     std::cout << "UNKNOWN TIME center=" << e.geometry().center()
    //               << " time =" << time
    //               << " starttime=" << starttime
    //               << " endtime=" << starttime+dt
    //               << std::endl;
    //   }
    
    // evaluate cell
    typename P0DGF::Traits::RangeType b(0.0);
    bathymmetry.evaluate(e,cellcenterlocal,b);
    typename P0DGF::Traits::RangeType uold(0.0);
    elevationold.evaluate(e,cellcenterlocal,uold);
    typename P0DGF::Traits::RangeType unew(0.0);
    elevationnew.evaluate(e,cellcenterlocal,unew);
    // std::cout << "STORAGE time=" << time << " starttime=" << starttime
    //           << " b=" << b[0] << " uold=" << uold[0] << " unew=" << unew[0]
    //           << " h=" << std::setprecision(9) << unew[0]-b[0]
    //           << " resh=" << residualheight << std::endl;

    // inactive cell
    if (unew[0]-b[0]<=residualheight)
      {
        // std::cout << "WET CELL detected time=" << time << " starttime=" << starttime << std::endl;
        if (std::abs(time-starttime)/dt<1e-6)
          return residualheight;
        else
          return residualheight;
      }

    // active cell
    auto w = (time-starttime)/dt;
    return std::max( (1.0-w)*(uold[0]-b[0])+w*(unew[0]-b[0]) ,0.0);
  }

  //! yes it is constant
  static constexpr bool permeabilityIsConstantPerCell()
  {
    return true;
  }

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return D;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename VDGF::Traits::RangeType velo(0.0);
    velocity.evaluate(e,x,velo);
    return velo;
  }

  //! sink term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! boundary condition type function
  /* return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet for Dirichlet boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann for flux boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow for outflow boundary conditions
   */
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& xlocal) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(xlocal);
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    return 1.0+0.8*std::sin(2*M_PI*xglobal[0]/10)*std::cos(2*M_PI*xglobal[1]/10);
  }

  //! flux boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! set time for subsequent evaluation
  void setTime (RF t) const
  {
    time = t;
    // std::cout << "modelT time set to " << time << std::endl;
  }

  void preStep (RF time_, RF dt_, int stages) const
  {
    starttime = time_;
    dt = dt_;
    // std::cout << "modelT preStep " << time << " " << dt << std::endl;
  }
};


//********************************************************************************
// driver function solving the problem
//********************************************************************************

template<typename GV>
void driverFV (const GV& gv, Dune::ParameterTree& ptree)
{
  // dimension and important types
  const int dim = GV::dimension;
  typedef typename GV::Grid::ctype DF;
  typedef double RF;                   // type for computations

  //=============================
  // set up the flow problem
  //=============================
  
  // Finite Element Map
  using FEM = Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim>;
  FEM fem(Dune::GeometryTypes::cube(dim));

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
  MODEL model;
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
  typedef typename GFS::template ConstraintsContainer<RF>::Type CC;
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
  int torder_flow = ptree.get("method.torder_flow",(int)1);
  Dune::PDELab::TimeSteppingParameterInterface<RF>*
    pmethod=&method1;
  if (torder_flow==1) pmethod = &method1;
  if (torder_flow==2) pmethod = &method2;
  if (torder_flow==3) pmethod = &method3;
  if (torder_flow<1||torder_flow>3)
    std::cout<<"torder not in [1,3]"<<std::endl;
  typedef Dune::PDELab::OneStepMethod<RF,IGO,PDESOLVER,Z,Z> OSM;
  OSM  osm(*pmethod,igo,pdesolver);
  osm.setVerbosityLevel(2);

  // velocity postprocessing
  Z znew(z);
  ZDGF znewdgf(gfs,znew);
  typedef NonlinearDiffusionFVVelocity<MODEL,GFS,Z> VeloDGF;
  auto residualheight = ptree.get<RF>("problem.residualheight");
  VeloDGF velodgf(model,gfs,znew,bathymmetry,residualheight);

  //=============================
  // set up the transport problem
  //=============================

  // model class
  using MODELT = ModelTransport<GV,ZDGF,VeloDGF>;
  MODELT modelT(ptree,bdgf,zdgf,znewdgf,velodgf);
  using GT=Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<MODELT>;
  GT gT(gv,modelT);
  using BT=Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<MODELT>;
  BT bT(modelT);
    
  // finite element map
  const int degree = 2;
  typedef Dune::PDELab::QkDGLocalFiniteElementMap<DF,RF,degree,dim,Dune::PDELab::QkDGBasisPolynomial::legendre> FEMT;
  const int nBlock = Dune::QkStuff::QkSize<degree,dim>::value;
  FEMT femT;

  // grid function space for transport problem
  typedef Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,nBlock> VBT;
  typedef Dune::PDELab::P0ParallelConstraints CONT;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEMT,CONT,VBT> GFST;
  GFST gfsT(gv,femT);
  gfsT.name("concentration");
  gfsT.ordering();
  std::cout << "concentration field has " << gfsT.globalSize() << " dofs" << std::endl;

  // coefficient vector
  using ZT = Dune::PDELab::Backend::Vector<GFST,RF>;
  ZT zT(gfsT);

  // make grid function for solution
  typedef Dune::PDELab::DiscreteGridFunction<GFST,ZT> DGFT;
  DGFT dgfT(gfsT,zT);

  // assemble constraints
  typedef typename GFST::template ConstraintsContainer<RF>::Type CCT;
  CCT ccT; ccT.clear();
  Dune::PDELab::constraints(bT,gfsT,ccT); // assemble constraints
  std::cout << "constrained dofs transport =" << ccT.size() << " of " << gfsT.globalSize() << std::endl;

  // initialize coefficient vector for transport with initial temperature field
  gT.setTime(time);
  Dune::PDELab::interpolate(gT,gfsT,zT);

  // assemblers for finite element problem
  typedef Dune::PDELab::ConvectionDiffusionDG<MODELT,FEMT> LOPT;
  LOPT lopT(modelT,Dune::PDELab::ConvectionDiffusionDGMethod::SIPG,Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn,3.0);
  typedef Dune::PDELab::GridOperator<GFST,GFST,LOPT,MBE,RF,RF,RF,CCT,CCT> GOST;
  GOST gosT(gfsT,ccT,gfsT,ccT,lopT,MBE(2*dim+1));
  //typedef Dune::PDELab::L2 MLOPT;
  typedef StorageL2<MODELT> MLOPT;
  MLOPT mlopT(modelT);
  typedef Dune::PDELab::GridOperator<GFST,GFST,MLOPT,MBE,RF,RF,RF,CCT,CCT> MGOST;
  MGOST mgosT(gfsT,ccT,gfsT,ccT,mlopT,MBE(2*dim+1));
  typedef Dune::PDELab::OneStepGridOperator<GOST,MGOST> IGOT;
  IGOT igoT(gosT,mgosT);
  //igoT.divideMassTermByDeltaT();
  igo.multiplySpatialTermByDeltaT();

  // linear solver backends
  //typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFST,CCT> LinearSolverT;
  typedef Dune::PDELab::ISTLBackend_OVLP_GMRES_ILU0<GFST,CCT> LinearSolverT;
  LinearSolverT lsT(gfsT,ccT,1000,2);

  // linear problem solver
  typedef Dune::PDELab::StationaryLinearProblemSolver<IGOT,LinearSolverT,ZT> PDESOLVERT;
  PDESOLVERT pdesolverT(igoT,lsT,1e-8);

  // implicit time-stepper
  int torder_transport = ptree.get("method.torder_transport",(int)1);
  Dune::PDELab::TimeSteppingParameterInterface<RF>* pmethodT;
  if (torder_transport==1) pmethodT = &method1;
  if (torder_transport==2) pmethodT = &method2;
  if (torder_transport==3) pmethodT = &method3;
  if (torder_transport<1||torder_transport>3)
    std::cout<<"torder not in [1,3]"<<std::endl;
  typedef Dune::PDELab::OneStepMethod<RF,IGOT,PDESOLVERT,ZT> OSMT;
  OSMT osmT(*pmethodT,igoT,pdesolverT);
  osmT.setVerbosityLevel(2);

  //=============================
  // visualization
  //=============================

  Z wetcells(gfs);
  ZDGF wetcellsdgf(gfs,wetcells);
  auto& indexset = gv.indexSet();
  for (const auto& e : elements(gv,Dune::Partitions::interior)) {
    auto i = indexset.index(e);
    auto bi = Dune::PDELab::Backend::native(bathymmetry)[i][0];
    auto zinew = Dune::PDELab::Backend::native(znew)[i][0];
    if (zinew-bi>residualheight)
      Dune::PDELab::Backend::native(wetcells)[i][0] = 1.0;
    else
      Dune::PDELab::Backend::native(wetcells)[i][0] = 0.0;
  }

  // prepare VTK writer and write first file
  std::string filename=ptree.get("output.filename","output");
  int subsampling=ptree.get("output.subsampling",(int)0);
  struct stat st;
  if( stat( filename.c_str(), &st ) != 0 )
    {
      int stat = 0;
      stat = mkdir(filename.c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
      if( stat != 0 && stat != -1)
        std::cout << "Error: Cannot create directory "
                  << filename << std::endl;
    }
  typedef Dune::SubsamplingVTKWriter<GV> VTKWRITER;
  VTKWRITER vtkwriter(gv,Dune::refinementLevels(subsampling));
  typedef Dune::VTKSequenceWriter<GV> VTKSEQUENCEWRITER;
  VTKSEQUENCEWRITER vtkSequenceWriter(std::make_shared<VTKWRITER>(vtkwriter),filename,filename,"");
  typedef Dune::PDELab::VTKGridFunctionAdapter<ZDGF> VTKF;
  vtkSequenceWriter.addCellData(std::shared_ptr<VTKF>(new VTKF(znewdgf,"solution")));
  vtkSequenceWriter.addCellData(std::shared_ptr<VTKF>(new VTKF(bdgf,"bathymmetry")));
  vtkSequenceWriter.addCellData(std::shared_ptr<VTKF>(new VTKF(wetcellsdgf,"wet cells")));
  typedef Dune::PDELab::VTKGridFunctionAdapter<VeloDGF> VeloVTKF;
  vtkSequenceWriter.addCellData(std::shared_ptr<VeloVTKF>(new VeloVTKF(velodgf,"velocity")));
  vtkSequenceWriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGFT> >(dgfT,"concentration"));
  vtkSequenceWriter.write(time,Dune::VTK::appendedraw);

  //=============================
  // time loop
  //=============================
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
      for (const auto& e : elements(gv,Dune::Partitions::interior)) {
        auto i = indexset.index(e);
        auto bi = Dune::PDELab::Backend::native(bathymmetry)[i][0];
        auto ziold = Dune::PDELab::Backend::native(z)[i][0];
        auto zinew = Dune::PDELab::Backend::native(znew)[i][0];
        max_change = std::max(max_change,std::abs(ziold-zinew));
        if (zinew-bi>residualheight)
          {
            Dune::PDELab::Backend::native(wetcells)[i][0] = 1.0;
          }
        else
          {
            Dune::PDELab::Backend::native(wetcells)[i][0] = 0.0;
            //Dune::PDELab::Backend::native(zT)[i] = 0.0;
          }
      }
      gv.comm().max(max_change);
      if (gv.comm().rank()==0)
        std::cout << "MAXCHANGE " << max_change << " " << time+dt << " " << dt << std::endl;

      // solve transport
      velodgf.update();  // note: velodgf depends on znew
      ZT zTnew(zT);
      osmT.apply(time,dt,zT,zTnew);
      
      // accept time step
      time+=dt;
      step++;
      z = znew;
      zT = zTnew;

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
    ptreeparser.readINITree("fourthexample.ini",ptree);
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

    driverFV(gv,ptree);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
