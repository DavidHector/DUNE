#ifndef Udune_ftworth_cdediscretizer_HH
#define Udune_ftworth_cdediscretizer_HH

#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/finiteelementmap/opbfem.hh>
#include<dune/pdelab/finiteelementmap/qkdg.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/constraints/p0.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/backend/istl.hh>

#include<dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include<dune/pdelab/localoperator/convectiondiffusiondg.hh>
#include<dune/pdelab/localoperator/convectiondiffusionccfv.hh>


// A base class for setting up a CDE problem with conforming finite elements
template<typename GV, typename Problem, typename FEM, int degree>
class CGConvectionDiffusionProblemBase
{
public:
  enum { dim = GV::dimension };
  using DF = typename GV::Grid::ctype;
  using RF = double;
  using BCFunction = Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem>;
  using G = Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem>;
  using CON = Dune::PDELab::ConformingDirichletConstraints;
  using VBE = Dune::PDELab::ISTL::VectorBackend<>;
  using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE>;
  using CC = typename GFS::template ConstraintsContainer<RF>::Type;
  using LOP = Dune::PDELab::ConvectionDiffusionFEM<Problem,FEM>;
  using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  using GO = Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,CC,CC>;
  using VEC = Dune::PDELab::Backend::Vector<GFS,RF>;
  using MAT = typename GO::Jacobian;
  using DGF = Dune::PDELab::DiscreteGridFunction<GFS,VEC>;

  const GV& gv;
  Problem& problem;
  std::shared_ptr<FEM> fem;
  std::shared_ptr<BCFunction> b;
  std::shared_ptr<G> g;
  std::shared_ptr<GFS> gfs;
  std::shared_ptr<CC> cc;
  std::shared_ptr<LOP> lop;
  std::shared_ptr<MBE> mbe;
  std::shared_ptr<GO> go;
  
  CGConvectionDiffusionProblemBase (const GV& gv_, Problem& problem_, const Dune::ParameterTree& ptree)
    : gv(gv_), problem(problem_)
  {
    // finite element map
    fem = std::shared_ptr<FEM>(new FEM(gv));

    // make parameter object and extract bctype and Dirichlet extension
    b = std::shared_ptr<BCFunction>(new BCFunction(gv,problem));
    g = std::shared_ptr<G>(new G(gv,problem));

    // Make grid function space
    gfs = std::shared_ptr<GFS>(new GFS(gv,*fem));
    gfs->name("Vh");

    // Assemble constraints
    cc = std::shared_ptr<CC>(new CC());
    Dune::PDELab::constraints(*b,*gfs,*cc); // assemble constraints
    std::cout << "constrained dofs=" << cc->size() << " of " << gfs->globalSize() << std::endl;

    // Make a local operator
    lop = std::shared_ptr<LOP>(new LOP(problem));

    // Make a global operator
    mbe = std::shared_ptr<MBE>(new MBE((int)pow((double)(1+2*degree),(double)dim)));
    go = std::shared_ptr<GO>(new GO(*gfs,*cc,*gfs,*cc,*lop,*mbe));
  }

};

template<typename GV, typename Problem, int degree>
class CGConvectionDiffusionProblemCube :
  public CGConvectionDiffusionProblemBase<GV,Problem,Dune::PDELab::QkLocalFiniteElementMap<GV,typename GV::Grid::ctype,double,degree>,degree>
{
  using FEM = Dune::PDELab::QkLocalFiniteElementMap<GV,typename GV::Grid::ctype,double,degree>;
  using Base = CGConvectionDiffusionProblemBase<GV,Problem,FEM,degree>;
public:
  CGConvectionDiffusionProblemCube (const GV& gv_, Problem& problem_, const Dune::ParameterTree& ptree)
    : Base(gv_,problem_,ptree)
  {}
};

template<typename GV, typename Problem, int degree>
class CGConvectionDiffusionProblemSimplex :
  public CGConvectionDiffusionProblemBase<GV,Problem,Dune::PDELab::PkLocalFiniteElementMap<GV,typename GV::Grid::ctype,double,degree>,degree>
{
  using FEM = Dune::PDELab::PkLocalFiniteElementMap<GV,typename GV::Grid::ctype,double,degree>;
  using Base = CGConvectionDiffusionProblemBase<GV,Problem,FEM,degree>;
public:
  CGConvectionDiffusionProblemSimplex (const GV& gv_, Problem& problem_, const Dune::ParameterTree& ptree)
    : Base(gv_,problem_,ptree)
  {}
};


// A base class for setting up a CDE problem with discontinuous Galerkin finite elements
template<typename GV, typename Problem, typename FEM, typename VBE, int degree>
class DGConvectionDiffusionProblemBase
{
public:
  enum { dim = GV::dimension };
  using DF = typename GV::Grid::ctype;
  using RF = double;
  using BCFunction = Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem>;
  using G = Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem>;
  using CON = Dune::PDELab::P0ParallelConstraints;
  using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE>;
  using CC = typename GFS::template ConstraintsContainer<RF>::Type;
  using LOP = Dune::PDELab::ConvectionDiffusionDG<Problem,FEM>;
  using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  using GO = Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,CC,CC>;
  using VEC = Dune::PDELab::Backend::Vector<GFS,RF>;
  using MAT = typename GO::Jacobian;
  using DGF = Dune::PDELab::DiscreteGridFunction<GFS,VEC>;

  const GV& gv;
  Problem& problem;
  std::shared_ptr<FEM> fem;
  std::shared_ptr<BCFunction> b;
  std::shared_ptr<G> g;
  std::shared_ptr<GFS> gfs;
  std::shared_ptr<CC> cc;
  std::shared_ptr<LOP> lop;
  std::shared_ptr<MBE> mbe;
  std::shared_ptr<GO> go;
  
  DGConvectionDiffusionProblemBase (const GV& gv_, Problem& problem_, const Dune::ParameterTree& ptree)
    : gv(gv_), problem(problem_)
  {
    // finite element map
    fem = std::shared_ptr<FEM>(new FEM());

    // make parameter object and extract bctype and Dirichlet extension
    b = std::shared_ptr<BCFunction>(new BCFunction(gv,problem));
    g = std::shared_ptr<G>(new G(gv,problem));

    // Make grid function space
    gfs = std::shared_ptr<GFS>(new GFS(gv,*fem));
    gfs->name("Vh");

    // Assemble constraints
    cc = std::shared_ptr<CC>(new CC());
    Dune::PDELab::constraints(*b,*gfs,*cc); // assemble constraints
    std::cout << "constrained dofs=" << cc->size() << " of " << gfs->globalSize() << std::endl;

    // Make a local operator
    lop = std::shared_ptr<LOP>(new LOP(problem,Dune::PDELab::ConvectionDiffusionDGMethod::SIPG,Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn,3.0));

    // Make a global operator
    mbe = std::shared_ptr<MBE>(new MBE(2*dim+1));
    go = std::shared_ptr<GO>(new GO(*gfs,*cc,*gfs,*cc,*lop,*mbe));
  }
};

template<typename GV, typename Problem, int degree>
class DGConvectionDiffusionProblemCube :
  public DGConvectionDiffusionProblemBase<GV,Problem,Dune::PDELab::QkDGLocalFiniteElementMap<typename GV::Grid::ctype,double,degree,GV::dimension,Dune::PDELab::QkDGBasisPolynomial::legendre>,Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,Dune::QkStuff::QkSize<degree,GV::dimension>::value>,degree>
{
  using FEM = Dune::PDELab::QkDGLocalFiniteElementMap<typename GV::Grid::ctype,double,degree,GV::dimension,Dune::PDELab::QkDGBasisPolynomial::legendre>;
  using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,Dune::QkStuff::QkSize<degree,GV::dimension>::value>;
  using Base = DGConvectionDiffusionProblemBase<GV,Problem,FEM,VBE,degree>;
public:
  DGConvectionDiffusionProblemCube (const GV& gv_, Problem& problem_, const Dune::ParameterTree& ptree)
    : Base(gv_,problem_,ptree)
  {}
};

template<typename GV, typename Problem, int degree>
class DGConvectionDiffusionProblemSimplex :
  public DGConvectionDiffusionProblemBase<GV,Problem,Dune::PDELab::OPBLocalFiniteElementMap<typename GV::Grid::ctype,double,degree,GV::dimension,Dune::GeometryType::simplex,Dune::GMPField<512> >,Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,Dune::PB::PkSize<degree,GV::dimension>::value>,degree>
{
  using FEM = Dune::PDELab::OPBLocalFiniteElementMap<typename GV::Grid::ctype,double,degree,GV::dimension,Dune::GeometryType::simplex,Dune::GMPField<512> >;
  using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,Dune::PB::PkSize<degree,GV::dimension>::value>;
  using Base = DGConvectionDiffusionProblemBase<GV,Problem,FEM,VBE,degree>;
public:
  DGConvectionDiffusionProblemSimplex (const GV& gv_, Problem& problem_, const Dune::ParameterTree& ptree)
    : Base(gv_,problem_,ptree)
  {}
};


// A base class for setting up a CDE problem with cell-centered finite volumes
template<typename GV, typename Problem, typename FEM>
class FVConvectionDiffusionProblemBase
{
public:
  enum { dim = GV::dimension };
  using DF = typename GV::Grid::ctype;
  using RF = double;
  using BCFunction = Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem>;
  using G = Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem>;
  using CON = Dune::PDELab::P0ParallelConstraints;
  using VBE = Dune::PDELab::ISTL::VectorBackend<>;
  using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE>;
  using CC = typename GFS::template ConstraintsContainer<RF>::Type;
  using LOP = Dune::PDELab::ConvectionDiffusionCCFV<Problem>;
  using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  using GO = Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,CC,CC>;
  using VEC = Dune::PDELab::Backend::Vector<GFS,RF>;
  using MAT = typename GO::Jacobian;
  using DGF = Dune::PDELab::DiscreteGridFunction<GFS,VEC>;

  const GV& gv;
  Problem& problem;
  std::shared_ptr<FEM> fem;
  std::shared_ptr<BCFunction> b;
  std::shared_ptr<G> g;
  std::shared_ptr<GFS> gfs;
  std::shared_ptr<CC> cc;
  std::shared_ptr<LOP> lop;
  std::shared_ptr<MBE> mbe;
  std::shared_ptr<GO> go;
  
  FVConvectionDiffusionProblemBase (const GV& gv_, Problem& problem_, const FEM& fem_, const Dune::ParameterTree& ptree)
    : gv(gv_), problem(problem_)
  {
    // finite element map
    fem = std::shared_ptr<FEM>(new FEM(fem_));

    // make parameter object and extract bctype and Dirichlet extension
    b = std::shared_ptr<BCFunction>(new BCFunction(gv,problem));
    g = std::shared_ptr<G>(new G(gv,problem));

    // Make grid function space
    gfs = std::shared_ptr<GFS>(new GFS(gv,*fem));
    gfs->name("Vh");

    // Assemble constraints
    cc = std::shared_ptr<CC>(new CC());
    Dune::PDELab::constraints(*b,*gfs,*cc); // assemble constraints
    std::cout << "constrained dofs=" << cc->size() << " of " << gfs->globalSize() << std::endl;

    // Make a local operator
    lop = std::shared_ptr<LOP>(new LOP(problem));

    // Make a global operator
    mbe = std::shared_ptr<MBE>(new MBE(2*dim+1));
    go = std::shared_ptr<GO>(new GO(*gfs,*cc,*gfs,*cc,*lop,*mbe));
  }
};

template<typename GV, typename Problem>
class FVConvectionDiffusionProblemCube :
  public FVConvectionDiffusionProblemBase<GV,Problem,Dune::PDELab::P0LocalFiniteElementMap<typename GV::Grid::ctype,double,GV::dimension>>
{
  using FEM = Dune::PDELab::P0LocalFiniteElementMap<typename GV::Grid::ctype,double,GV::dimension>;
  using Base = FVConvectionDiffusionProblemBase<GV,Problem,FEM>;
public:
  FVConvectionDiffusionProblemCube (const GV& gv_, Problem& problem_, const Dune::ParameterTree& ptree)
    : Base(gv_,problem_,FEM(Dune::GeometryTypes::cube(GV::dimension)),ptree)
  {}
};


#endif // Udune_ftworth_HH
