// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/*************************************/
/* This example solves surface flow in the vicinity of Heidelberg
   with 30m resolution and constant rainfall. Uses structured mesh
   and finite volume method.
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
#include<dune/pdelab/solver/newton.hh>

// dune-vtk includes
#include <dune/vtk/writers/vtkimagedatawriter.hh>
#include <dune/vtk/writers/vtkrectilineargridwriter.hh>
#include <dune/vtk/writers/vtkstructuredgridwriter.hh>
#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>
#include <dune/vtk/datacollectors/yaspdatacollector.hh>
#include <dune/vtk/pvdwriter.hh>

// include stuff from dune-hydro
#include<dune/hydro/nonlineardiffusionfv.hh>
#include<dune/hydro/nonlineardiffusionfv_velocity.hh>
#include<dune/hydro/geotiffreader.hh>
#include<dune/hydro/netcdfreader.hh>
#include<dune/hydro/rasterdataset.hh>

// a function for outlier correction
template<typename T>
RasterDataSet<T> outlier (RasterDataSet<T>& image_in, int diameter, int threshold)
{
  // prepare result
  RasterDataSet<T> image_out(image_in);
  std::vector<T>& out = image_out.data();
  diameter = std::max(1,diameter);

  // analyse size etc of input image
  auto& in=image_in.data(); // access to raw image data
  int N=in.size(); // number of pixels in image
  int m=image_in.sizeLong(); // pixels per row
  int n=image_in.sizeLat(); // pixels per column

  // compute average value of pixels
  for (int j=0; j<n; j++)
    for (int i=0; i<m; i++)
      {
        std::vector<T> values;
        for (int jj=std::max(j-diameter,0); jj<std::min(j+diameter+1,n); jj++)
          for (int ii=std::max(i-diameter,0); ii<std::min(i+diameter+1,m); ii++)
            values.push_back(in[jj*m+ii]);
        std::sort(values.begin(),values.end());
        T median;
        if (values.size()%2==0)
          median = (values[values.size()/2]+values[values.size()/2-1])/2;
        else
          median = values[values.size()/2];
        if (std::abs(median-in[j*m+i])>=threshold)
          {
            std::cout << "i=" << i << " j=" << j << " value=" << in[j*m+i] << " median=" << median << std::endl;
            out[j*m+i] = median;
          }
        else
          out[j*m+i] = in[j*m+i];
      }
  return image_out;
}

//********************************************************************************
/** \brief Model class provides parameters to a shallow flow model

    The PDE is:
     \partial_t(s(x) u) - \nabla \cdot (k(x,u,b) phi(u,b,\nabla u) \nabla u) = f(x,t) in \Omega
                                                                           u = u0(x)  at t = 0
                                                                           u = g(x,t) on \Gamma_D
                              -(k(x,u,b) phi(u,b,\nabla u) \nabla u)\cdot\nu = j(x,t) on \Gamma_N

   with
     u(x,t)             position of free surface of water body [m]
     b(x)               position of bottom [m]
     s(x)               storage term [-]
     k(x,u,b)           permeability
     phi(u,b,\nabla u)  nonlinearity
     f                  source/sink term
     g                  Dirichlet boundary condition
     j                  flux boundary condition
     u0                 initial value

   \tparam Number : a number type
 */
//********************************************************************************
template<typename Number, typename DGF>
class Model
{
  enum {dim=2};
  Dune::ParameterTree ptree;
  Dune::FieldVector<Number,dim> L;
  Dune::FieldVector<Number,dim> H;
  std::array<int,dim> N;
  Number precipitationrate;
  Number initialheight;
  RasterDataSet<short> elevation_raster;
  std::shared_ptr<DGF> psourcedgf;
  Number eps,A,B;

public:
  //! Constructor gets parameter tree to read ini file parameters
  Model (Dune::ParameterTree ptree_)
    : ptree(ptree_), eps(1e-2)
  {
    // The Baden-Württemberg data
    GeoTIFFImage<GInt16> hs_con_n45e005("hs_con_n45e005.tif","/Users/peterbastian/Data/BWat90m/",1,2);
    GeoTIFFImage<GInt16> hs_con_n45e010("hs_con_n45e010.tif","/Users/peterbastian/Data/BWat90m/",1,2);

    // The Baden-Württemberg data
    double dx=std::abs(hs_con_n45e005.dLong());
    double dy=std::abs(hs_con_n45e005.dLat());


    double ox = ptree.get("problem.ox",(double) 7.0);
    double oy = ptree.get("problem.oy",(double)46.0);
    N[0] = ptree.get("problem.NX",(int)4800);
    N[1] = ptree.get("problem.NY",(int)4800);
    H[0] = 90.0;
    H[1] = 90.0;
    L[0] = N[0]*H[0];
    L[1] = N[1]*H[1];
    elevation_raster = RasterDataSet<short>(ox+0.5*dx,oy+0.5*dy,dx,dy,N[0],N[1],0.0,1);
    elevation_raster.paste(hs_con_n45e005);
    elevation_raster.paste(hs_con_n45e010);
    elevation_raster = outlier(elevation_raster,2,2000);

    // read precipitation rate
    precipitationrate = ptree.get("problem.precipitationrate",(Number)0.0);
    precipitationrate = precipitationrate/86.4e6; // convert from liter per day per m^2 to meter per second
    initialheight = ptree.get("problem.initialheight",(Number)0.0);

    // regularization parameter
    A = 2.5/(2.0*std::sqrt(eps));
    B = 0.5/(2.0*pow(eps,2.5));
    std::cout << "REGULARIZATION gradient eps=" << eps << " A=" << A << " B=" << B << std::endl;
  }

  void set_source (std::shared_ptr<DGF> p)
  {
    psourcedgf = p;
  }
  
  const Dune::FieldVector<Number,dim>& length () const
  {
    return L;
  }

  const std::array<int,dim>& cells () const
  {
    return N;
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
  template<typename E, typename X>
  /**
     \param e      : codim 0 entity
     \param xlocal : local coordinate
   */
  Number bathymmetry (const E& e, const X& xlocal) const
  {
    auto x = e.geometry().center();
    int cellx = std::floor((x[0]-0.5*H[0])/H[0]);
    int celly = std::floor((x[1]-0.5*H[1])/H[1]);
    return elevation_raster(cellx,celly);
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
  Number phi_1 (U b_inside, U b_outside, U u_inside, U u_outside, VN vn) const
  {
    // upwind evaluation of height
    Number u_upwind = (vn>=0) ? u_inside : u_outside;
    
    // evaluation of (u-b)
    Number bmax=std::max(b_inside,b_outside);
    Number h_upwind = std::max(u_upwind-bmax,0.0);
    
    return h_upwind;
  }

  template<typename U, typename VN>
  Number dphi_du_inside_1 (U b_inside, U b_outside, U u_inside, U u_outside, VN vn) const
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
  Number dphi_du_outside_1 (U b_inside, U b_outside, U u_inside, U u_outside, VN vn) const
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
  Number dphi_dvn_1 (U b_inside, U b_outside, U u_inside, U u_outside, VN vn) const
  {
    return 0.0;
  }

  //----------------
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
    
    // regularization of gradient and evaluation
    Number sign = (vn>=0) ? 1.0 : -1.0;
    if (std::abs(vn)>=eps)
      return -h_upwind*0.5/(std::abs(vn)*std::sqrt(std::abs(vn)))*sign;
    else
      return -h_upwind*B*2.0*vn;
  }
  //----------------

  
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
    using Y = Dune::FieldVector<Number,1>;
    Y result;
    psourcedgf->evaluate(e,x,result);
    return result[0]*precipitationrate;
  }

  //! boundary condition type function (true = Dirichlet)
  /**
     \param i intersection (with boundary)
     \param x local coordinate on intersection
   */
  template<typename I, typename X>
  bool b (const I& i, const X& x) const
  {
    return true;
  }

  //! Dirichlet extension
  /** Provides Dirichlet boundary values as well as initial value!
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number g (const E& e, const X& xlocal) const
  {
    auto height = bathymmetry(e,xlocal);
    auto x = e.geometry().global(xlocal);
    Number eps=1e-4;
    if (x[0]>eps && x[0]<L[0]-eps && x[1]>eps && x[1]<L[1]-eps)
      return height+initialheight;
    else
      return height;
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

void driverFVNewNewton (Dune::ParameterTree& ptree)
{
  // dimension and important types
  const int dim = 2;
  using RF = double;               // type for computations

  using Grid = Dune::YaspGrid<dim>;
  using DF = Grid::ctype;
  using GV = Grid::LeafGridView;

  using FEM = Dune::PDELab::P0LocalFiniteElementMap<RF,RF,dim>;
  using CON = Dune::PDELab::P0ParallelConstraints;
  using VBE = Dune::PDELab::ISTL::VectorBackend<>;
  using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE>;

  using Z = Dune::PDELab::Backend::Vector<GFS,RF>;

  using ZDGF = Dune::PDELab::DiscreteGridFunction<GFS,Z>;
  using VTKF = Dune::PDELab::VTKGridFunctionAdapter<ZDGF>;

  using MODEL = Model<RF,ZDGF>;

  // make time
  RF time = 0.0;

  // make model
  MODEL model(ptree);
  model.setTime(time);

  // make YaspGrid
  std::shared_ptr<Grid> gridp = std::make_shared<Grid>(model.length(),model.cells(),std::bitset<dim>(0ULL),1);
  GV gv=gridp->leafGridView();

  // Make grid function space
  FEM fem(Dune::GeometryTypes::cube(dim));
  CON con;
  GFS gfs(gv,fem,con);
  gfs.name("Vh");

  // A coefficient vector
  Z z(gfs); // initial value

  // Make a grid function out of it
  ZDGF zdgf(gfs,z);

  // make user functions and set initial time
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
  LS ls (gfs,400,1);
  auto params = ls.parameters();
  params.setCoarsenTarget(80000);
  ls.setParameters(params);
  
  // using LS = Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,CC>;
  // LS ls(gfs,cc,500,5,2);

  // solve nonlinear problem
  typedef Dune::PDELab::NewtonMethod<IGO,LS> PDESOLVER;
  PDESOLVER pdesolver(igo,ls);
  pdesolver.setParameters(ptree.sub("newton"));
  if (gv.comm().rank()!=0)
    pdesolver.setVerbosityLevel(0);

  pdesolver.setReassembleThreshold(0.0);
  pdesolver.setVerbosityLevel(2);
  pdesolver.setAbsoluteLimit(1e-7);
  pdesolver.setReduction(1e-6);
  pdesolver.setMinLinearReduction(1e-4);
  pdesolver.setUseMaxNorm(true);
  //pdesolver.setMaxIterations(10);
  //pdesolver.setLineSearchMaxIterations(8);
  //pdesolver.setLineSearchStrategy("hackbuschReuskenAcceptBest");
  //pdesolver.setLineSearchStrategy("hackbuschReusken");

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
  using VeloDGF = NonlinearDiffusionFVVelocity<MODEL,GFS,Z>;
  using VeloVTKF = Dune::PDELab::VTKGridFunctionAdapter<VeloDGF>;
  VeloDGF velodgf(model,gfs,znew,bathymmetry);

  // residuals
  Z rspace(gfs); rspace = 0.0;  ZDGF rspacedgf(gfs,rspace);
  Z rtime(gfs); rtime = 0.0;    ZDGF rtimedgf(gfs,rtime);
  Z rtotal(gfs); rtotal = 0.0;  ZDGF rtotaldgf(gfs,rtotal);
  Z zdiff(gfs); zdiff = 0.0;    ZDGF zdiffdgf(gfs,zdiff);

  // source term
  auto initialheight = ptree.get<double>("problem.initialheight");
  Z source(gfs);
  auto psourcedgf = std::make_shared<ZDGF>(gfs,source);
  for (int i=0; i<Dune::PDELab::Backend::native(z).N(); ++i)
    Dune::PDELab::Backend::native(source)[i][0] = std::max(0.0,1.0-(Dune::PDELab::Backend::native(z)[i][0]-Dune::PDELab::Backend::native(bathymmetry)[i][0])/initialheight);
  model.set_source(psourcedgf);

  // prepare VTK writer and write first file
  std::string filename=ptree.get("output.filename","output");
  // struct stat st;
  // if( stat( filename.c_str(), &st ) != 0 )
  //   {
  //     int stat = 0;
  //     stat = mkdir(filename.c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
  //     if( stat != 0 && stat != -1)
  //       std::cout << "Error: Cannot create directory "
  //                 << filename << std::endl;
  //   }

  // new writer
  using Writer = Dune::VtkImageDataWriter<GV>;
  Dune::PvdWriter<Writer> pvdWriter(gv, Dune::Vtk::FormatTypes::COMPRESSED, Dune::Vtk::DataTypes::FLOAT32);
  pvdWriter.addCellData(std::make_shared<VTKF>(znewdgf,"solution"));
  pvdWriter.addCellData(std::make_shared<VTKF>(bdgf,"bathymmetry"));
  pvdWriter.addCellData(std::make_shared<VTKF>(*psourcedgf,"source"));
  pvdWriter.addCellData(std::make_shared<VeloVTKF>(velodgf,"velocity"));
  std::string fullfilename = filename + ".vti";
  pvdWriter.writeTimestep(time,fullfilename);
  
  // time loop
  RF T = ptree.get("problem.T",(RF)1.0);
  RF dt = ptree.get("method.dt",(RF)0.1);
  RF timestepmax =  ptree.get("method.dtmax",dt);
  int every = ptree.get("output.every",(int)1);
  int step=0;
  int increased_step = 0;
  int decreased_step = 0;
  double eps=1e-6;
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
      catch (Dune::PDELab::LineSearchError) {
        dt *= factor;
        if (dt<eps) throw;
        decreased_step = step;
        if (gv.comm().rank()==0)
          std::cout << "STEP CONTROL: " << " decreased in step " << step+1 << " newdt: " << dt << std::endl;
        continue;
      }
      catch (Dune::PDELab::TerminateError) {
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

      // analyze residual
      rspace = 0.0;
      go0.residual(znew,rspace);
      Dune::PDELab::set_constrained_dofs(cc,0.0,rspace);
      rtime=0.0;
      zdiff=znew;
      zdiff -= z;
      go1.residual(zdiff,rtime);
      Dune::PDELab::set_constrained_dofs(cc,0.0,rtime);
      //rtime *= 1.0/dt; // this seems to be there already!
      rtotal = rtime;
      rtotal += rspace;

      auto norm2_time = std::sqrt(gv.comm().sum(rtime.two_norm2()));
      auto norm2_space = std::sqrt(gv.comm().sum(rspace.two_norm2()));
      auto norm2_total = std::sqrt(gv.comm().sum(rtotal.two_norm2()));
      if (gv.comm().rank()==0)
        {
          std::cout << "||rtime ||_2=" << norm2_time  << std::endl;
          std::cout << "||rspace||_2=" << norm2_space  << std::endl;
          std::cout << "||rtotal||_2=" << norm2_total  << std::endl;
        }
      auto norm8_time = gv.comm().max(rtime.infinity_norm());
      auto norm8_space = gv.comm().max(rspace.infinity_norm());
      auto norm8_total = gv.comm().max(rtotal.infinity_norm());
      if (gv.comm().rank()==0)
        {
          std::cout << "||rtime ||_8=" << norm8_time  << std::endl;
          std::cout << "||rspace||_8=" << norm8_space  << std::endl;
          std::cout << "||rtotal||_8=" << norm8_total  << std::endl;
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

      // update source term
      for (int i=0; i<Dune::PDELab::Backend::native(z).N(); ++i)
        Dune::PDELab::Backend::native(source)[i][0] = std::max(0.0,1.0-(Dune::PDELab::Backend::native(z)[i][0]-Dune::PDELab::Backend::native(bathymmetry)[i][0])/initialheight);

      // determine maximum velicity value in any edge
      auto maxv = maxvelocity(velodgf);
      if (gv.comm().rank()==0)
        std::cout << "MAXVELOCITY " << maxv << std::endl;
      
      // output to VTK file
      if (step%every==0)
        {
          pvdWriter.writeTimestep(time,fullfilename);
          // vtkSequenceWriter.write(time,Dune::VTK::appendedraw);
          if (gv.comm().rank()==0)
            std::cout << "WRITING VTK OUTPUT " << step << " " << time << " " << dt << std::endl;
        }

      // increase timestep
      if (step-decreased_step>24  && step-increased_step>24)
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
    ptreeparser.readINITree("bw90mpre.ini",ptree);
    ptreeparser.readOptions(argc,argv,ptree);

    driverFVNewNewton(ptree);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
