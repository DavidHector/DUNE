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
#include"spe10.hh"

#define DEGREE 1
#define USECG 1
//#define USEDG 1
//#define USEFV 1
#define USEMLGENEO 1
//#define USEAMG 1
//#define USEILU 1

/*! Adapter that extracts diagonal of permeability tensor from parameter class

  \tparam T  model of ConvectionDiffusionParameterInterface
*/
template<typename T>
class DiagonalPermeabilityAdapter
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
									   typename T::Traits::RangeFieldType,
									   T::Traits::dimDomain,Dune::FieldVector<typename T::Traits::RangeFieldType,T::Traits::dimDomain> >
					  ,DiagonalPermeabilityAdapter<T> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
					   typename T::Traits::RangeFieldType,
					   T::Traits::dimDomain,Dune::FieldVector<typename T::Traits::RangeFieldType,T::Traits::dimDomain> > Traits;

  //! constructor
  DiagonalPermeabilityAdapter (const typename Traits::GridViewType& g_, T& t_)
    : g(g_), t(t_)
  {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
			const typename Traits::DomainType& x,
			typename Traits::RangeType& y) const
  {
    for (int i=0; i<T::Traits::dimDomain; i++)
      y[i] = log10(t.A(e,x)[i][i]);
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return g;
  }

  inline void setTime(double time_)
  {
    t.setTime(time_);
  }

private:
  const typename Traits::GridViewType& g;
  T& t;
};

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
    ptreeparser.readINITree("spe10.ini",ptree);
    ptreeparser.readOptions(argc,argv,ptree);

    //===================
    // grid generation and partitioning
    //===================
    
    // make a YaspGrid
    using Grid = Dune::YaspGrid<3>;  
    Dune::FieldVector<double,3> upper_right;
    upper_right[0] = ptree.get("grid.LX",(double)1.0);
    upper_right[1] = ptree.get("grid.LY",(double)1.0);
    upper_right[2] = ptree.get("grid.LZ",(double)1.0);
    std::array<int,3> cells;
    cells[0] = ptree.get("grid.NX",(int)10);
    cells[1] = ptree.get("grid.NY",(int)10);
    cells[2] = ptree.get("grid.NZ",(int)10);
    std::shared_ptr<Grid> pgrid = std::shared_ptr<Grid>(new Grid(upper_right,cells));

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

    size_t overlap = ptree.get("geneo.overlap",(size_t)1);
    std::string extensionmethod = ptree.get("geneo.extensionmethod","vertex");
    size_t drop_small_overlap = ptree.get("geneo.drop_small_overlap",(size_t)1);
    bool coordinate_partitioning = ptree.get("geneo.coordinate_partitioning",(bool)false);
    std::vector<size_t> parts = ptree.get<std::vector<size_t>>("geneo.parts");

    // make domain decomposition on fine grid
    std::shared_ptr<DomainDecomposition> pdd;
    if (coordinate_partitioning)
      pdd = std::make_shared<DomainDecomposition>(gv,subdomains[0],parts,overlap,extensionmethod,drop_small_overlap,false);
    else
      pdd = std::make_shared<DomainDecomposition>(gv,Dune::MPIHelper::getLocalCommunicator(),subdomains[0],overlap,extensionmethod,drop_small_overlap,false);

    //===================
    // set up a discretization scheme
    //===================

    // solve a Poisson problem
    using DF = Grid::ctype;
    using RF = double;

    // set up the solver class
    using Problem = SPE10Problem<GV,RF>;
    Problem problem(ptree);
#ifdef USECG
    using Discretization = CGConvectionDiffusionProblemCube<GV,Problem,DEGREE>;
    bool hasSkeleton = false;
#endif
#ifdef USEDG
    using Discretization = DGConvectionDiffusionProblemCube<GV,Problem,DEGREE>;
    bool hasSkeleton = true;
#endif
#ifdef USEFV
    using Discretization = FVConvectionDiffusionProblemCube<GV,Problem>;
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
    // choose random initial guess
    for (size_t i=0; i<Dune::PDELab::Backend::native(x).N(); i++)
      for (size_t j=0; j<Dune::PDELab::Backend::native(x)[i].N(); j++)
        Dune::PDELab::Backend::native(x)[i][j] = (rand() % 10000)/10000.0;
    set_constrained_dofs(*disc.cc,0.0,x); // dirichlet dofs are zero, others are 1, so we can use it as a mask!
    MAT A(*disc.go); // this makes the sparsity pattern
    A = 0.0;
    ISTLM& istlA = Dune::PDELab::Backend::native(A);
    std::shared_ptr<ISTLM> pistlA = stackobject_to_shared_ptr(Dune::PDELab::Backend::native(istlA));

    // residual and correction
    VEC r(*disc.gfs);
    r = 0.0;
    disc.go->residual(x,r);
    ISTLV& istlr = Dune::PDELab::Backend::native(r);
    
    //===================
    // build subdomain information & test matrices in subsets
    //==================
    
    auto pdddm = std::shared_ptr<DomainDecompositionDOFMapper>(new DomainDecompositionDOFMapper(*disc.gfs,pdd));

    std::cout << "create sparse matrices for snippets and subdomains" << std::endl;

    auto pvolume_snippet_matrices = std::shared_ptr<std::vector<std::shared_ptr<ISTLM>>>( new std::vector<std::shared_ptr<ISTLM>>(set_up_local_matrices(istlA,pdddm->volumesnippetnumber_local_to_global,pdddm->volumesnippetnumber_global_to_local)) );
    for (auto matptr : *pvolume_snippet_matrices) *matptr = 0.0; // clear subdomain matrices

    std::shared_ptr<std::vector<std::shared_ptr<ISTLM>>> pskeleton_snippet_matrices; // the empty pointer
    if (hasSkeleton)
      {
        pskeleton_snippet_matrices = std::shared_ptr<std::vector<std::shared_ptr<ISTLM>>>( new std::vector<std::shared_ptr<ISTLM>>(set_up_local_matrices(istlA,pdddm->skeletonsnippetnumber_local_to_global,pdddm->skeletonsnippetnumber_global_to_local)) );
        for (auto matptr : *pskeleton_snippet_matrices) *matptr = 0.0; // clear subdomain matrices
      }

    //===================
    // matrix assembly
    //==================

    std::cout << "assembling subdomain matrices" << std::endl;
    // now we may wrap the local operator
    using WLOP = SnippetAssemblerWrapper<typename Discretization::LOP,ISTLM,typename Discretization::GFS>;
    WLOP wlop(*disc.gfs,*disc.lop,pdddm,pvolume_snippet_matrices,pskeleton_snippet_matrices);

    // and set up a new grid operator
    using WGOP = Dune::PDELab::GridOperator<typename Discretization::GFS,typename Discretization::GFS,WLOP,
                                            typename Discretization::MBE,typename Discretization::RF,typename Discretization::RF,typename Discretization::RF,
                                            typename Discretization::CC,typename Discretization::CC>;
    WGOP wgop(*disc.gfs,*disc.cc,*disc.gfs,*disc.cc,wlop,*disc.mbe);
    Dune::Timer timer;
    timer.reset();
    wgop.jacobian(x,A); // assembles the global stiffness matrix and the snippets
                        // we have global dirichlet constraints in the global matrix A but *not* in the snippets!
    auto assemble_time = timer.elapsed();
    std::cout << "ASSEMBLE TIME=" << assemble_time << std::endl;
    
    // global vector with dirichlet boundary conditions
    // contains 0.0 for dof on Dirichlet boundary, otherwise 1.0
    VEC dirichlet_mask(*disc.gfs);
    dirichlet_mask = 1.0;
    set_constrained_dofs(*disc.cc,0.0,dirichlet_mask); // dirichlet dofs are zero, others are 1, so we can use it as a mask!
    std::shared_ptr<ISTLV> pistldirichlet_mask = stackobject_to_shared_ptr(Dune::PDELab::Backend::native(dirichlet_mask));


    //===================
    // now set up solver and solve
    //==================
    
    size_t coarse_overlap = ptree.get("geneo.coarse_overlap",(size_t)0);
    std::string pum = ptree.get("geneo.pum","standard");
    std::string fineGEVPrhs = ptree.get("geneo.fineGEVPrhs","pu");
    std::string coarseGEVPrhs = ptree.get("geneo.coarseGEVPrhs","pu");
    std::string coarseeigensolver = ptree.get("geneo.coarseeigensolver","eigen");
    std::string cycle = ptree.get("geneo.cycle","additive");
    size_t n_eigenvectors_fine_computed = ptree.get("geneo.n_eigenvectors_fine_computed",(size_t)10);
    size_t n_eigenvectors_fine_used = ptree.get("geneo.n_eigenvectors_fine_used",(size_t)5);
    size_t n_eigenvectors_coarse_computed = ptree.get("geneo.n_eigenvectors_coarse_computed",(size_t)10);
    size_t n_eigenvectors_coarse_used = ptree.get("geneo.n_eigenvectors_coarse_used",(size_t)5);
    size_t mlgeneo_verbose = ptree.get("geneo.verbose",(size_t)0);
    double abs_zero_ker = ptree.get("geneo.abs_zero_ker",(double)1e-10);
    double regularization_ker = ptree.get("geneo.regularization_ker",(double)0.0);
    double eigenvalue_fine_threshold = ptree.get("geneo.eigenvalue_fine_threshold",(double)-1.0);
    double eigenvalue_coarse_threshold = ptree.get("geneo.eigenvalue_coarse_threshold",(double)-1.0);
    double arpack_tolerance = ptree.get("geneo.arpack_tolerance",(double)0.0);
    bool merge_disconnected = ptree.get("geneo.merge_disconnected",(bool)false);

    using Operator = Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV>;
    Operator op(istlA);
    std::cout << "creating multilevel geneo preconditioner" << std::endl;
    timer.reset();
#ifdef USEMLGENEO
    MultiLevelGenEO<ISTLM,ISTLV> prec(Dune::MPIHelper::getLocalCommunicator(),
                                      pistlA,
                                      pistldirichlet_mask,
                                      pvolume_snippet_matrices,
                                      pskeleton_snippet_matrices,
                                      pdddm,
                                      subdomains,
                                      pum,
                                      fineGEVPrhs,
                                      coarseGEVPrhs,
                                      coarse_overlap,
                                      coarseeigensolver,
                                      arpack_tolerance,
                                      n_eigenvectors_fine_computed,n_eigenvectors_fine_used,
                                      n_eigenvectors_coarse_computed,n_eigenvectors_coarse_used,
                                      eigenvalue_fine_threshold,
                                      eigenvalue_coarse_threshold,
                                      abs_zero_ker,
                                      regularization_ker,
                                      merge_disconnected,
                                      cycle,
                                      mlgeneo_verbose);
#endif
#ifdef USEAMG    
    using Smoother = Dune::SeqSSOR<ISTLM,ISTLV,ISTLV>;
    using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments ;
    using AMG = Dune::Amg::AMG<Operator,ISTLV,Smoother> ;
    using Parameters = Dune::Amg::Parameters;
    using Criterion = Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<ISTLM,Dune::Amg::FirstDiagonal>> ;
    Parameters params(15,2000);
    params.setDefaultValuesIsotropic(dim);
    params.setDebugLevel(1);
    SmootherArgs smootherArgs;
    smootherArgs.iterations = 1;
    smootherArgs.relaxationFactor = 1;
    Criterion criterion(params);
    AMG prec(op,criterion,smootherArgs);
#endif
#ifdef USEILU
    size_t ilun = ptree.get("ilu.n",(size_t)0);
    double iluw = ptree.get("ilu.w",(double)1.0);
    Dune::SeqILU<ISTLM,ISTLV,ISTLV> prec(istlA,ilun,iluw,false);
#endif
    auto setup_time = timer.elapsed();
    std::cout << "SETUP TIME=" << setup_time << std::endl;
    size_t itmax = ptree.get("solver.itmax",(size_t)100);
    std::string solvertype = ptree.get("solver.type","gmres");
    Dune::CGSolver<ISTLV> cgsolver(op,prec,1e-8,itmax,2);
    Dune::RestartedGMResSolver<ISTLV> gmressolver(op,prec,1e-8,50,2500,2);
    Dune::InverseOperatorResult stat;
    ISTLV istlv(istlx);
    istlv = (double)0.0;
    if (solvertype=="gmres")
      gmressolver.apply(istlv,istlr,stat);
    else
      cgsolver.apply(istlv,istlr,stat);      
    istlx -= istlv;

    //===================
    // visualization
    //==================

    // generate some vectors for visualization
#ifdef USEMLGENEO
    std::vector<std::vector<unsigned>> k(subdomains.size());
    for (auto& v : k) v.resize(indexset.size(0));
    size_t k_0 = 0;
    unsigned count_elements=0;
    for (size_t i=0; i<indexset.size(0); ++i)
      {
        for (size_t l=0; l<subdomains.size(); l++)
          k[l][i] = prec.get_overlappingsubdomains(l)[i].size();
        k_0 = std::max(k_0,prec.get_overlappingsubdomains(0)[i].size());
        count_elements += prec.get_overlappingsubdomains(0)[i].size();
      }
    std::cout << k_0 << " k_0" << std::endl;
    std::cout << indexset.size(0) << " elements, " << count_elements << " in all subdomains, ratio=" << ((double)count_elements)/indexset.size(0) << std::endl;
#endif

    bool save_output = ptree.get("output.save",(bool)true);
    if (!save_output) return 0; // then we are done

    // visualize result
    std::string outfile = ptree.get("output.filename","output");
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::nonconforming);
    using DGF = typename Discretization::DGF;
    DGF dgfx(*disc.gfs,x);
    using VTKF = Dune::PDELab::VTKGridFunctionAdapter<DGF>;
    vtkwriter.addVertexData(std::shared_ptr<VTKF>(new VTKF(dgfx,"fesol")));
    typedef DiagonalPermeabilityAdapter<Problem> PermDGF;
    PermDGF permdgf(gv,problem);
    typedef Dune::PDELab::VTKGridFunctionAdapter<PermDGF> PermVTKDGF;
    vtkwriter.addCellData(std::make_shared<PermVTKDGF>(permdgf,"K"));
#ifdef USEMLGENEO
    vtkwriter.addCellData(pdd->partition,"partition 0");
    vtkwriter.addCellData(pdd->vizpartition,"vizpartition 0");
    vtkwriter.addCellData(prec.get_partition(1),"partition 1");
    vtkwriter.addCellData(prec.get_vizpartition(1),"vizpartition 1");
    for (size_t l=0; l<subdomains.size(); l++)
      {
        std::ostringstream s1;
        s1 << "volumesnippetnumber" << l ;
        vtkwriter.addCellData(prec.get_element_to_volumesnippetnumber(l),s1.str());
      }
    for (size_t l=0; l<subdomains.size(); l++)
      {
        std::ostringstream s1;
        s1 << "k" << l ;
        vtkwriter.addCellData(k[l],s1.str());
      }
    bool vizbasisvectors = ptree.get("geneo.vizbasisvectors",false);
    int n_bv=0;
    std::vector<std::vector<size_t>> viz_vectors;
    if (vizbasisvectors)
      n_bv = viz_vectors.size();
    std::vector<std::shared_ptr<VEC>> basisvectors(n_bv); // empty vector
    std::vector<std::shared_ptr<DGF>> basisvectordgfs(n_bv); // the discrete grid functions
    std::vector<std::string> basisvectornames(n_bv); // the names
    if (vizbasisvectors)
      for (int i=0; i<n_bv; i++)
        {
          int l = viz_vectors[i][0];
          int p = viz_vectors[i][1];
          int k = viz_vectors[i][2];
          basisvectors[i] = std::shared_ptr<VEC>( new VEC(*disc.gfs,0.0) ); // add a vector
          Dune::PDELab::Backend::native(*basisvectors[i]) = prec.prolongated_basis_vector(l,p,k); // fill value
          basisvectordgfs[i] = std::shared_ptr<DGF>( new DGF(*disc.gfs,*basisvectors[i]) );
          std::ostringstream s1;
          s1 << "basisvector_" << l << "_" << p << "_" << k;
          basisvectornames[i] = s1.str();
        }
    for (int i=0; i<basisvectors.size(); i++)
      {
        std::cout << "add output for " << basisvectornames[i] << std::endl;
        vtkwriter.addVertexData(std::shared_ptr<VTKF>(new VTKF(*basisvectordgfs[i],basisvectornames[i])));
      }
#endif
    vtkwriter.write(outfile,Dune::VTK::appendedraw);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
