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
#include<dune/pdelab/solver/newton.hh>

// include stuff from dune-hydro
#include<dune/hydro/nonlineardiffusionfv.hh>
#include<dune/hydro/nonlineardiffusionfv_velocity.hh>

//GDAL includes. We require that GDAL is found
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()

//********************************************************************************
// GDAL stuff. to be put into a seperate file later ...
//********************************************************************************

/** list all tiles needed for a certain region
    The arguments are given in latitude and longitude:
    x_lower   longitude of most west point, west of Greenwich is negative, east is positive
    y_lower   latitude of most south point, south of equator is negative, north of equator is positive
    x_upper   longitude of most east point, west of Greenwich is negative, east is positive
    y_upper   latitude of most north point, south of equator is negative, north of equator is positive
 */
void listTiles (double x_lower, double y_lower, double x_upper, double y_upper)
{
  // tiles have the format [n|s]YY_[e,w]XXX_1arc_v3.tif
  int xint_lower = floor(x_lower);
  int yint_lower = floor(y_lower);
  int xint_upper = ceil(x_upper);
  int yint_upper = ceil(y_upper);
  for (int y=yint_lower; y<yint_upper; y++)
    for (int x=xint_lower; x<xint_upper; x++)
      {
        std::stringstream latstring;
        if (y>=0) latstring << "n" << y; else latstring << "s" << std::abs(y);      
        std::stringstream longstring;
        if (x>=0) longstring << "e" << x; else longstring << "w" << std::abs(x);
        std::stringstream tilename;
        tilename << latstring.str() << "_" << longstring.str() << "_1arc_v3.tif"; 
        std::cout << "need tile " << x << "," << y << " to " << x+1 << "," << y+1 << "  " << tilename.str() << std::endl;
      }
}

/** check if all tiles are available
 */
bool checkTiles ( double x_lower, double y_lower, double x_upper, double y_upper, std::string path)
{
}

/** A class that reads a GeoTIFFImage
    
    - holds the data in a buffer and lets you access the pixels
    - provides size of the image
    - provides positions of pixel centers in lat-long
    - Note: images are y-downwards!
 */
class GeoTIFFImage {
public:
  GeoTIFFImage (std::string filename, std::string path="")
  {
    GDALDataset  *poDataset;
    poDataset = (GDALDataset *) GDALOpen( filename.c_str(), GA_ReadOnly );
    if ( poDataset == NULL )
      {
        std::cout << "GDALOpen() failed" << std::endl;
        // Need to throw an exception here
      }

    // read size
    xsize = poDataset->GetRasterXSize();
    ysize = poDataset->GetRasterYSize();

    // read transformation to lat long coordinates
    double GT[6];
    if ( poDataset->GetGeoTransform( GT ) != CE_None )
      {
        std::cout << "GetGeoTransform() failed" << std::endl;
        // Need to throw an exception here
      }
    oX = GT[0];  oY = GT[3];
    AXX = GT[1]; AXY = GT[2];
    AYX = GT[4]; AYY = GT[5];
    
    // read file contents
    GDALRasterBand* band = poDataset->GetRasterBand(1);
    int nXBlockSize, nYBlockSize;
    band->GetBlockSize( &nXBlockSize, &nYBlockSize );
    int nXBlocks = (band->GetXSize() + nXBlockSize - 1) / nXBlockSize;
    int nYBlocks = (band->GetYSize() + nYBlockSize - 1) / nYBlockSize;
    p = new GInt16[xsize*ysize]; // allocate image
    GInt16 *pabyData = new GInt16[nXBlocks*nYBlocks]; // allocate transfer buffer

    // loop over blocks
    for ( int iYBlock = 0; iYBlock < nYBlocks; iYBlock++ )
      for ( int iXBlock = 0; iXBlock < nXBlocks; iXBlock++ )
        {
          auto x=band->ReadBlock( iXBlock, iYBlock, pabyData );
          int nXValid, nYValid;
          band->GetActualBlockSize(iXBlock, iYBlock, &nXValid, &nYValid); // infer valid part of buffer
          for (int iy=0; iy<nYValid; iy++)
            { // loop over valid pixels in buffer
              auto iImage = ((iYBlock*nYBlockSize+iy)*xsize)+(iXBlock*nXBlockSize); // first pixel of line in image
              auto iBuffer = iy*nXBlockSize; // first pixel of line in buffer
              for (int ix=0; ix<nXValid; ix++) p[iImage++] = pabyData[iBuffer++];
            }
        }

    delete[] pabyData; // delete buffer
  }

  GeoTIFFImage (const GeoTIFFImage& image)
  {
    // deep copy
    xsize = image.xsize;
    ysize = image.ysize;
    
    oX = image.oX;
    oY = image.oY;
    AXX = image.AXX;
    AXY = image.AXY;
    AYX = image.AYX;
    AYY = image.AYY;
    
    p = new GInt16[xsize*ysize]; // allocate image
    for (int i=0; i<xsize*ysize; i++) p[i] = image.p[i];
  }
  
  // destructor
  ~GeoTIFFImage ()
  {
    delete[] p;
  }
  
  // number of pixels in x (longitude)
  int sizeX () const
  {
    return xsize;
  }

  // number of pixels in y (latitude)
  int sizeY () const
  {
    return ysize;
  }

  // longitude of pixel center
  double longitude (int x, int y) const
  {
    return oX + AXX*(x+0.5) + AXY*(y+0.5);
  }

  // latitude of pixel center
  double latitude (int x, int y) const
  {
    return oY + AYX*(x+0.5) + AYY*(y+0.5);
  }

  // read pixel 
  GInt16 getPixel (int x, int y) const
  {
    return p[y*xsize+x];
  }

  // read pixel 
  void putPixel (int x, int y, GInt16 value)
  {
    p[y*xsize+x] = value;
  }

  // simplifed access operator with round brackets
  const GInt16& operator() (int x, int y) const
  {
    return p[y*xsize+x];
  }
  GInt16& operator() (int x, int y)
  {
    return p[y*xsize+x];
  }

private:
  int xsize, ysize; // number of pixels in x and y
  double oX,oY;     // origin of image
  double AXX, AXY, AYX, AYY; // transformation matrix
  GInt16 *p; // image buffer
};

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
template<typename Number>
class Model
{
  Number eps, A,B;
  Dune::ParameterTree ptree;
  Dune::FieldVector<Number,2> L;
  std::array<int,2> N;
  Number precipitationrate;
  std::vector<std::vector<float>> heightdata;
  int CellsX[5];
  int CellsY[5];
  int level,ox,oy;

public:
  //! Constructor gets parameter tree to read ini file parameters
  Model (Dune::ParameterTree ptree_)
    : eps(1e-3), ptree(ptree_), heightdata(5)
  {
    // read global domain and mesh size
    L[0] = ptree.get("grid.structured.LX",(Number)1.0);
    L[1] = ptree.get("grid.structured.LY",(Number)1.0);
    N[0] = ptree.get("grid.structured.NX",(unsigned int) 10);
    N[1] = ptree.get("grid.structured.NY",(unsigned int) 10);
    precipitationrate = ptree.get("problem.precipitationrate",(Number)0.0);

    // load an image
    GeoTIFFImage image("n49_e008_1arc_v3.tif");
    if (image.sizeX()-1!=3600 || image.sizeY()-1!=3600)
      {
        std::cout << "only works with 3600^2 cells" << std::endl;
        return;
      }

    // allocate a hierarchy of representations
    int Lmax = 4;
    CellsX[0] = CellsY[0] = 225;
    for (int l=1; l<=Lmax; l++) CellsX[l] = CellsY[l] = 2*CellsX[l-1];
    for (int l=Lmax; l>=0; l--) heightdata[l].resize(CellsX[l]*CellsY[l]);

    // compute coarser representations of height
    for (int y=0; y<CellsY[Lmax]; y++)
      for (int x=0; x<CellsX[Lmax]; x++)
        heightdata[Lmax][y*CellsX[Lmax]+x] = image(x,CellsY[Lmax]-1-y);
    for (int l=Lmax; l>0; l--)
      {
        for (int i=0; i< heightdata[l-1].size(); i++)
          heightdata[l-1][i] = 1.0;
        for (int y=0; y<CellsY[l]; y++)
          for (int x=0; x<CellsX[l]; x++)
            heightdata[l-1][(y/2)*CellsX[l-1]+(x/2)] *= heightdata[l][y*CellsX[l]+x];
        for (int i=0; i< heightdata[l-1].size(); i++)
          heightdata[l-1][i] = std::pow(heightdata[l-1][i],0.25);
      }

    // level to use in height data hierarchy
    level = ptree.get("bathymmetry.level",(int) 0);
    ox = ptree.get("bathymmetry.ox",(int) 0);
    oy = ptree.get("bathymmetry.oy",(int) 0);
    std::cout << "using level=" << level << " ox=" << ox << " oy=" << oy << std::endl; 

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
  template<typename E, typename X>
  /**
     \param e      : codim 0 entity
     \param xlocal : local coordinate
   */
  Number bathymmetry (const E& e, const X& xlocal) const
  {
    auto x = e.geometry().center();
    // determine cell number in x and y
    int cellx = std::floor(x[0]/(L[0]/N[0])) + ox;
    int celly = std::floor(x[1]/(L[1]/N[1])) + oy;

    int pixel = celly*CellsX[level]+cellx;
    //std::cout << cellx << " " << celly << " | " << CellsX[level] << " " << CellsY[level] << " | " << pixel << std::endl;
    
    return heightdata[level][pixel];
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
    return precipitationrate;
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
    auto global = i.geometry().global(x);
    if (global[0]<1e-6 && std::abs(global[1]-50.0)<5.0)
      return true;
  }

  //! Dirichlet extension
  /** Provides Dirichlet boundary values as well as initial value!
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number g (const E& e, const X& x) const
  {
    auto b = bathymmetry(e,x);
    return b+1e-3;
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

  // residuals
  Z rspace(gfs); rspace = 0.0;  ZDGF rspacedgf(gfs,rspace);
  Z rtime(gfs); rtime = 0.0;    ZDGF rtimedgf(gfs,rtime);
  Z rtotal(gfs); rtotal = 0.0;  ZDGF rtotaldgf(gfs,rtotal);
  Z zdiff(gfs); zdiff = 0.0;    ZDGF zdiffdgf(gfs,zdiff);

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
  vtkSequenceWriter.addCellData(std::shared_ptr<VTKF>(new VTKF(rspacedgf,"rspace")));
  vtkSequenceWriter.addCellData(std::shared_ptr<VTKF>(new VTKF(rtimedgf,"rtime")));
  vtkSequenceWriter.addCellData(std::shared_ptr<VTKF>(new VTKF(rtotaldgf,"rtotal")));
  vtkSequenceWriter.addCellData(std::shared_ptr<VTKF>(new VTKF(zdiffdgf,"zdiff")));
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


template<typename GV, typename FEM>
void driverFVNewNewton (const GV& gv, const FEM& fem, Dune::ParameterTree& ptree)
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
  typedef Dune::PDELab::NewtonMethod<IGO,LS> PDESOLVER;
  PDESOLVER pdesolver(igo,ls);
  pdesolver.setReassembleThreshold(0.0);
  pdesolver.setVerbosityLevel(2);
  pdesolver.setAbsoluteLimit(1e-10);
  pdesolver.setReduction(1e-8);
  pdesolver.setMinLinearReduction(1e-4);
  pdesolver.setUseMaxNorm(true);
  //pdesolver.setMaxIterations(15);
  //pdesolver.setLineSearchMaxIterations(10);
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
  typedef NonlinearDiffusionFVVelocity<MODEL,GFS,Z> VeloDGF;
  VeloDGF velodgf(model,gfs,znew,bathymmetry);

  // residuals
  Z rspace(gfs); rspace = 0.0;  ZDGF rspacedgf(gfs,rspace);
  Z rtime(gfs); rtime = 0.0;    ZDGF rtimedgf(gfs,rtime);
  Z rtotal(gfs); rtotal = 0.0;  ZDGF rtotaldgf(gfs,rtotal);
  Z zdiff(gfs); zdiff = 0.0;    ZDGF zdiffdgf(gfs,zdiff);

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
  vtkSequenceWriter.addCellData(std::shared_ptr<VTKF>(new VTKF(rspacedgf,"rspace")));
  vtkSequenceWriter.addCellData(std::shared_ptr<VTKF>(new VTKF(rtimedgf,"rtime")));
  vtkSequenceWriter.addCellData(std::shared_ptr<VTKF>(new VTKF(rtotaldgf,"rtotal")));
  vtkSequenceWriter.addCellData(std::shared_ptr<VTKF>(new VTKF(zdiffdgf,"zdiff")));
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
    ptreeparser.readINITree("thirdexample.ini",ptree);
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
    std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(new Grid(L,N,std::bitset<dim>(0ULL),1));

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

    //typedef Dune::PDELab::QkLocalFiniteElementMap<GV,double,double,1> FEM;
    // typedef Dune::PDELab::PkLocalFiniteElementMap<GV,double,double,1> FEM;
    // FEM fem(gv);
    // driverFEM(gv,fem,ptree);

    typedef Dune::PDELab::P0LocalFiniteElementMap<double,double,dim> FEM;
    FEM fem(Dune::GeometryTypes::cube(dim));
    driverFVNewNewton(gv,fem,ptree);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
