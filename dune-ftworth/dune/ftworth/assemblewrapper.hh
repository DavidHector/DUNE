#ifndef Udune_ftworth_assemblewrapper_HH
#define Udune_ftworth_assemblewrapper_HH

#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/callswitch.hh>

#include"subdomainutilities.hh"

/** \brief A local operator that accumulates contributions to the subdomain matrices
 *
 * Assumptions:
 * - Problem is linear to begin with
 * - wrapper should only be used for calling jacobian on grid operator
 * - no thought went into the apply methods yet
 *
 * Notes:
 * - Still need to think about how to handle Dirichlet boundary conditions on the domain boundary
 *
 * \tparam[in] LocalOperator Type of the local operator that gets wrapped
 */
template <typename LocalOperator, typename ISTLM, typename FLSDI>
class SubdomainAssemblerWrapper
  : public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
  // copy flags of wrapped operator
  static constexpr bool doPatternVolume = LocalOperator::doPatternVolume;
  static constexpr bool doPatternVolumePostSkeleton = LocalOperator::doPatternVolumePostSkeleton;
  static constexpr bool doPatternSkeleton = LocalOperator::doPatternSkeleton;
  static constexpr bool doPatternBoundary = LocalOperator::doPatternBoundary;

  static constexpr bool doAlphaVolume = LocalOperator::doAlphaVolume;
  static constexpr bool doAlphaVolumePostSkeleton = true; //LocalOperator::doAlphaVolumePostSkeleton;
  static constexpr bool doAlphaSkeleton = LocalOperator::doAlphaSkeleton;
  static constexpr bool doAlphaBoundary = LocalOperator::doAlphaBoundary;

  // lambda methods
  static constexpr bool doLambdaVolume = LocalOperator::doLambdaVolume;
  static constexpr bool doLambdaVolumePostSkeleton = LocalOperator::doLambdaVolumePostSkeleton;
  static constexpr bool doLambdaSkeleton = LocalOperator::doLambdaSkeleton;
  static constexpr bool doLambdaBoundary = LocalOperator::doLambdaBoundary;

  // If the underlying lop is linear, this one will be linear too
  static constexpr bool isLinear = LocalOperator::isLinear;

  // copy face handling
  static constexpr bool doSkeletonTwoSided = LocalOperator::doSkeletonTwoSided;

  /** \brief Construct new instance of class
   *
   * \param[in] _localOperator Wrapped local operator instance
   */
  SubdomainAssemblerWrapper(const LocalOperator& localOperator_,
			    const std::vector<std::set<size_t>>& overlappingsubdomains_,
			    FLSDI& flsdi_,
			    std::shared_ptr<std::vector<std::shared_ptr<ISTLM>>> matrices_)
    : localOperator(localOperator_), overlappingsubdomains(overlappingsubdomains_), flsdi(flsdi_), matrices(matrices_)
  {}

  /** \brief Copy constructor */
  SubdomainAssemblerWrapper(const SubdomainAssemblerWrapper& other)
    : localOperator(other.localOperator), overlappingsubdomains(other.overlappingsubdomains), flsdi(other.flsdi), matrices(other.matrices)
  {}

  //////////
  // pattern methods
  //////////
      
  // define sparsity pattern on element
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume (const LFSU& lfsu, const LFSV& lfsv,
                       LocalPattern& pattern) const
  {
    Dune::PDELab::LocalOperatorApply::patternVolume(localOperator,lfsu,lfsv,pattern);
  }

  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume_post_skeleton (const LFSU& lfsu, const LFSV& lfsv,
				     LocalPattern& pattern) const
  {
    Dune::PDELab::LocalOperatorApply::patternVolumePostSkeleton(localOperator,lfsu,lfsv,pattern);
  }

  // define sparsity pattern connecting self and neighbor dofs
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_skeleton (const LFSU& lfsu_s, const LFSV& lfsv_s, const LFSU& lfsu_n, const LFSV& lfsv_n,
                         LocalPattern& pattern_sn,
                         LocalPattern& pattern_ns) const
  {
    Dune::PDELab::LocalOperatorApply::patternSkeleton(localOperator,lfsu_s,lfsv_s,lfsu_n,lfsv_n,pattern_sn,pattern_ns);
  }


  // define sparsity pattern connecting dofs on boundary elements
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_boundary(const LFSU& lfsu_s, const LFSV& lfsv_s,
                        LocalPattern& pattern_ss) const
  {
    Dune::PDELab::LocalOperatorApply::patternBoundary(localOperator,lfsu_s,lfsv_s,pattern_ss);
  }

  //////////
  // alpha methods
  //////////

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    Dune::PDELab::LocalOperatorApply::alphaVolume(localOperator,eg,lfsu,x,lfsv,r);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume_post_skeleton (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    Dune::PDELab::LocalOperatorApply::alphaVolumePostSkeleton(localOperator,eg,lfsu,x,lfsv,r);
  }

  // skeleton integral depending on test and ansatz functions
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_skeleton (const IG& ig,
                       const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                       const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                       R& r_s, R& r_n) const
  {
    Dune::PDELab::LocalOperatorApply::alphaSkeleton(localOperator,ig,lfsu_s,x_s,lfsv_s,lfsu_n,x_n,lfsv_n,r_s,r_n);
  }

  // boundary integral depending on test and ansatz functions
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig,
                       const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                       R& r_s) const
  {
    Dune::PDELab::LocalOperatorApply::alphaBoundary(localOperator,ig,lfsu_s,x_s,lfsv_s,r_s);
  }

  //////////
  // lambda methods
  //////////

  // volume integral depending only on test functions
  template<typename EG, typename LFSV, typename R>
  void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
  {
    Dune::PDELab::LocalOperatorApply::lambdaVolume(localOperator,eg,lfsv,r);
  }

  template<typename EG, typename LFSV, typename R>
  void lambda_volume_post_skeleton (const EG& eg, const LFSV& lfsv, R& r) const
  {
    Dune::PDELab::LocalOperatorApply::lambdaVolumePostSkeleton(localOperator,eg,lfsv,r);
  }
  template<typename IG, typename LFSV, typename R>
  void lambda_skeleton(const IG& ig,
		       const LFSV& lfsv_s, const LFSV& lfsv_n,
		       R& r_s, R& r_n) const
  {
    Dune::PDELab::LocalOperatorApply::lambdaSkeleton(localOperator,ig,lfsv_s,lfsv_n,r_s,r_n);
  }
  template<typename IG, typename LFSV, typename R>
  void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
  {
    Dune::PDELab::LocalOperatorApply::lambdaBoundary(localOperator,ig,lfsv,r);
  }

  //////////
  // jacobian methods
  //////////
      
  template<typename EG, typename LFSU, typename X, typename LFSV, typename MAT>
  void jacobian_volume (const EG& eg,
                        const LFSU& lfsu,
                        const X& x,
                        const LFSV& lfsv,
                        MAT& mat) const
  {
    // std::cout << "jacobian_volume on " << flsdi.indexset.index(eg.entity()) << std::endl;
    
    // call local operator
    Dune::PDELab::LocalOperatorApply::jacobianVolume(localOperator,eg,lfsu,x,lfsv,mat);

    // the scatter is done in jacobian_volume_post_skeleton
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename MAT>
  void jacobian_volume_post_skeleton (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, MAT& mat) const
  {
    //std::cout << "jacobian_volume_post_skeleton on " << flsdi.indexset.index(eg.entity()) << std::endl;

    // call local operator
    Dune::PDELab::LocalOperatorApply::jacobianVolumePostSkeleton(localOperator,eg,lfsu,x,lfsv,mat);

    // load information about the indices
    flsdi.lfs.bind(eg.entity());
    flsdi.lfscache.update();
    std::vector<unsigned> localdata(flsdi.lfs.size());
    flsdi.view.bind(flsdi.lfscache);
    flsdi.view.read(localdata);
    flsdi.view.unbind();

    // get underlying container for local stiffness matrix
    const auto& container = mat.container();

    // scatter local stiffness matrix into the subdomain matrices
    // this is the scatter for all contributions to the diagonal block
    for (auto p : overlappingsubdomains[flsdi.indexset.index(eg.entity())])
      {
	//std::cout << "element " << flsdi.indexset.index(eg.entity()) << " in subdomain " << p << std::endl;
	// now we know that subdomain p contains this element
	for (unsigned i = 0; i<flsdi.lfs.size(); ++i)
	  {
	    auto globalblocki = localdata[i]/flsdi.blocksize;
	    auto compi = localdata[i]%flsdi.blocksize;
	    auto subdomainblocki = (*flsdi.global_to_local)[globalblocki].find(p)->second; 
	    //std::cout << "   globalblocki, compi, subdomainblocki : " << globalblocki << " " << compi << " " << subdomainblocki << std::endl;
	    for (unsigned j = 0; j<flsdi.lfs.size(); ++j)
	      {
		auto globalblockj = localdata[j]/flsdi.blocksize;
		auto compj = localdata[j]%flsdi.blocksize;
		auto subdomainblockj = (*flsdi.global_to_local)[globalblockj].find(p)->second;
		//std::cout << "      globalblockj, compj, subdomainblockj : " << globalblockj << " " << compj << " " << subdomainblockj << std::endl;
		(*(*matrices)[p])[subdomainblocki][subdomainblockj][compi][compj] += container(lfsv,i,lfsu,j); 
	      }
	  }
      }
  }

  template<typename IG, typename LFSU, typename X, typename LFSV, typename MAT>
  void jacobian_skeleton (const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                          MAT& mat_ss, MAT& mat_sn,
                          MAT& mat_ns, MAT& mat_nn) const
  {
    //std::cout << "jacobian_skeleton on " << flsdi.indexset.index(ig.inside()) << "," << flsdi.indexset.index(ig.outside()) << std::endl;

    // call local operator
    auto container_ss_before = mat_ss.container(); // copy container, so we may compute what has changed
    Dune::PDELab::LocalOperatorApply::jacobianSkeleton(localOperator,ig,lfsu_s,x_s,lfsv_s,lfsu_n,x_n,lfsv_n,mat_ss,mat_sn,mat_ns,mat_nn);

    // load information about the indices
    flsdi.lfs.bind(ig.inside());
    flsdi.lfscache.update();
    std::vector<unsigned> localdata_s(flsdi.lfs.size());
    flsdi.view.bind(flsdi.lfscache);
    flsdi.view.read(localdata_s);
    flsdi.view.unbind();

    flsdi.lfs.bind(ig.outside());
    flsdi.lfscache.update();
    std::vector<unsigned> localdata_n(flsdi.lfs.size());
    flsdi.view.bind(flsdi.lfscache);
    flsdi.view.read(localdata_n);
    flsdi.view.unbind();

    // get underlying container for local stiffness matrix
    const auto& container_ss = mat_ss.container();
    const auto& container_sn = mat_sn.container();
    const auto& container_ns = mat_ns.container();
    const auto& container_nn = mat_nn.container();
    
    // scatter local stiffness matrix contributions into the subdomain matrices
    // only if a subdomain contains *both* elements (and then all four matrices must be scattered
    for (auto p : overlappingsubdomains[flsdi.indexset.index(ig.inside())])
      if (overlappingsubdomains[flsdi.indexset.index(ig.outside())].count(p)>0)
	{
	  // now we know that subdomain p contains both elements
	  // we also assume that degrees of freedom for both elements as well as ansatz and test space are the same
	  for (unsigned i = 0; i<flsdi.lfs.size(); ++i)
	    {
	      auto globalblocki_s = localdata_s[i]/flsdi.blocksize;
	      auto globalblocki_n = localdata_n[i]/flsdi.blocksize;
	      auto compi_s = localdata_s[i]%flsdi.blocksize; // same because block size is same
	      auto compi_n = localdata_n[i]%flsdi.blocksize; // same because block size is same
	      auto subdomainblocki_s = (*flsdi.global_to_local)[globalblocki_s].find(p)->second; 
	      auto subdomainblocki_n = (*flsdi.global_to_local)[globalblocki_n].find(p)->second; 
	      for (unsigned j = 0; j<flsdi.lfs.size(); ++j)
		{
		  auto globalblockj_s = localdata_s[j]/flsdi.blocksize;
		  auto globalblockj_n = localdata_n[j]/flsdi.blocksize;
		  auto compj_s = localdata_s[j]%flsdi.blocksize;
		  auto compj_n = localdata_n[j]%flsdi.blocksize;
		  auto subdomainblockj_s = (*flsdi.global_to_local)[globalblockj_s][p];
		  auto subdomainblockj_n = (*flsdi.global_to_local)[globalblockj_n][p];
		  (*(*matrices)[p])[subdomainblocki_s][subdomainblockj_n][compi_s][compj_n] += container_sn(lfsu_s,i,lfsu_n,j); 
		  (*(*matrices)[p])[subdomainblocki_n][subdomainblockj_s][compi_n][compj_s] += container_ns(lfsu_n,i,lfsu_s,j); 
	          (*(*matrices)[p])[subdomainblocki_n][subdomainblockj_n][compi_n][compj_n] += container_nn(lfsu_n,i,lfsu_n,j); 
		}
	    }
	}
      else
      	{
      	  // now here only subdomain p contains the inside element but not the outside element
	  // SO WE MAY NOT USE global_to_local in subdomain p for the outside element!
      	  // we also assume that degrees of freedom for both elements as well as ansatz and test space are the same
      	  for (unsigned i = 0; i<flsdi.lfs.size(); ++i)
      	    {
      	      auto globalblocki_s = localdata_s[i]/flsdi.blocksize;
      	      auto compi_s = localdata_s[i]%flsdi.blocksize; // same because block size is same
      	      auto subdomainblocki_s = (*flsdi.global_to_local)[globalblocki_s].find(p)->second; 
      	      for (unsigned j = 0; j<flsdi.lfs.size(); ++j)
      		{
      		  auto globalblockj_s = localdata_s[j]/flsdi.blocksize;
      		  auto compj_s = localdata_s[j]%flsdi.blocksize;
      		  auto subdomainblockj_s = (*flsdi.global_to_local)[globalblockj_s][p];
      		  (*(*matrices)[p])[subdomainblocki_s][subdomainblockj_s][compi_s][compj_s] -= container_ss(lfsu_s,i,lfsu_s,j)-container_ss_before(lfsu_s,i,lfsu_s,j); 
      		}
      	    }    
      	}
  }

  template<typename IG, typename LFSU, typename X, typename LFSV, typename MAT>
  void jacobian_boundary (const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          MAT& mat_ss) const
  {
    // std::cout << "jacobian_boundary on " << flsdi.indexset.index(ig.inside()) << std::endl;

    // call local operator
    Dune::PDELab::LocalOperatorApply::jacobianBoundary(localOperator,ig,lfsu_s,x_s,lfsv_s,mat_ss);

    // scatter is done in post skeleton!
  }

  //////////
  // jacobian apply methods
  //////////


private:
  /** \brief Wrapped original local operator */
  const LocalOperator& localOperator;
  const std::vector<std::set<size_t>>& overlappingsubdomains;
  FLSDI& flsdi;
  std::shared_ptr<std::vector<std::shared_ptr<ISTLM>>> matrices;
};



/** \brief A local operator that accumulates contributions to the snippet matrices
 *
 * Assumptions:
 * - Problem is linear to begin with
 * - wrapper should only be used for calling jacobian on grid operator
 * - no thought went into the apply methods yet
 *
 * Notes:
 * - Still need to think about how to handle Dirichlet boundary conditions on the domain boundary
 *
 * \tparam[in] LocalOperator Type of the local operator that gets wrapped
 */
template <typename LocalOperator, typename ISTLM, typename GFS>
class SnippetAssemblerWrapper
  : public Dune::PDELab::LocalOperatorDefaultFlags
{
  // some types
  using GV = typename GFS::Traits::GridView;
  using IS = typename GV::IndexSet;
  using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
  using LFSCache = Dune::PDELab::LFSIndexCache<LFS>;
  using RF = unsigned int;
  using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
  using ISTLV = Dune::PDELab::Backend::Native<Z>;
  using VIEW = typename Z::template LocalView<LFSCache>;

  // local data
  const LocalOperator& localOperator; // Wrapped original local operator
  std::shared_ptr<DomainDecompositionDOFMapper> pdddm;
  std::shared_ptr<std::vector<std::shared_ptr<ISTLM>>> pvolume_snippet_matrices;
  std::shared_ptr<std::vector<std::shared_ptr<ISTLM>>> pskeleton_snippet_matrices;
  const GV& gv; // from the gfs
  const IS& indexset; // from the gv
  mutable LFS lfs; // a local function space for the gfs
  mutable LFSCache lfscache; // lfs cache for accessing dofs of an element
  Z dofnumber; // one-dimensional flat index for each dof
  mutable VIEW view; // local part of the vector associated with an element
  size_t blocksize; // assume that vector consists of blocks of this fixed size
  DomainDecompositionDOFMapper& dddm;
  DomainDecomposition& dd;

public:
  // copy flags of wrapped operator
  static constexpr bool doPatternVolume = LocalOperator::doPatternVolume;
  static constexpr bool doPatternVolumePostSkeleton = LocalOperator::doPatternVolumePostSkeleton;
  static constexpr bool doPatternSkeleton = LocalOperator::doPatternSkeleton;
  static constexpr bool doPatternBoundary = LocalOperator::doPatternBoundary;

  static constexpr bool doAlphaVolume = LocalOperator::doAlphaVolume;
  static constexpr bool doAlphaVolumePostSkeleton = LocalOperator::doAlphaVolumePostSkeleton;
  static constexpr bool doAlphaSkeleton = LocalOperator::doAlphaSkeleton;
  static constexpr bool doAlphaBoundary = LocalOperator::doAlphaBoundary;

  // lambda methods
  static constexpr bool doLambdaVolume = LocalOperator::doLambdaVolume;
  static constexpr bool doLambdaVolumePostSkeleton = LocalOperator::doLambdaVolumePostSkeleton;
  static constexpr bool doLambdaSkeleton = LocalOperator::doLambdaSkeleton;
  static constexpr bool doLambdaBoundary = LocalOperator::doLambdaBoundary;

  // If the underlying lop is linear, this one will be linear too
  static constexpr bool isLinear = LocalOperator::isLinear;

  // copy face handling
  static constexpr bool doSkeletonTwoSided = LocalOperator::doSkeletonTwoSided;

  /** \brief Construct new instance of class
   *
   * \param[in] _localOperator Wrapped local operator instance
   */
  SnippetAssemblerWrapper(const GFS& gfs,
			  const LocalOperator& localOperator_,
			  std::shared_ptr<DomainDecompositionDOFMapper> pdddm_,
			  std::shared_ptr<std::vector<std::shared_ptr<ISTLM>>> pvolume_snippet_matrices_,
			  std::shared_ptr<std::vector<std::shared_ptr<ISTLM>>> pskeleton_snippet_matrices_
			  )
    : localOperator(localOperator_), pdddm(pdddm_), pvolume_snippet_matrices(pvolume_snippet_matrices_), pskeleton_snippet_matrices(pskeleton_snippet_matrices_),
      gv(gfs.gridView()), indexset(gv.indexSet()), lfs(gfs), lfscache(lfs), dofnumber(gfs), view(dofnumber), dddm(*pdddm_), dd(*pdddm->pdd)
  {
    // first we fill a vector with the global dof numbers
    // we assume that our ISTL vector is a BlockVector<FieldVector<>>
    RF dof_counter = 0; 
    ISTLV& istldofnumber = Dune::PDELab::Backend::native(dofnumber); // a vector storing the global dof number in each dof
    for (int i=0; i<istldofnumber.size(); i++)
      for (int j=0; j<istldofnumber[i].size(); j++)
        istldofnumber[i][j] = dof_counter++;

    // determine the block size
    blocksize = istldofnumber[0].size();
  }

  //////////
  // pattern methods
  //////////
      
  // define sparsity pattern on element
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume (const LFSU& lfsu, const LFSV& lfsv,
                       LocalPattern& pattern) const
  {
    Dune::PDELab::LocalOperatorApply::patternVolume(localOperator,lfsu,lfsv,pattern);
  }

  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume_post_skeleton (const LFSU& lfsu, const LFSV& lfsv,
				     LocalPattern& pattern) const
  {
    Dune::PDELab::LocalOperatorApply::patternVolumePostSkeleton(localOperator,lfsu,lfsv,pattern);
  }

  // define sparsity pattern connecting self and neighbor dofs
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_skeleton (const LFSU& lfsu_s, const LFSV& lfsv_s, const LFSU& lfsu_n, const LFSV& lfsv_n,
                         LocalPattern& pattern_sn,
                         LocalPattern& pattern_ns) const
  {
    Dune::PDELab::LocalOperatorApply::patternSkeleton(localOperator,lfsu_s,lfsv_s,lfsu_n,lfsv_n,pattern_sn,pattern_ns);
  }


  // define sparsity pattern connecting dofs on boundary elements
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_boundary(const LFSU& lfsu_s, const LFSV& lfsv_s,
                        LocalPattern& pattern_ss) const
  {
    Dune::PDELab::LocalOperatorApply::patternBoundary(localOperator,lfsu_s,lfsv_s,pattern_ss);
  }

  //////////
  // alpha methods
  //////////

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    Dune::PDELab::LocalOperatorApply::alphaVolume(localOperator,eg,lfsu,x,lfsv,r);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume_post_skeleton (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    Dune::PDELab::LocalOperatorApply::alphaVolumePostSkeleton(localOperator,eg,lfsu,x,lfsv,r);
  }

  // skeleton integral depending on test and ansatz functions
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_skeleton (const IG& ig,
                       const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                       const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                       R& r_s, R& r_n) const
  {
    Dune::PDELab::LocalOperatorApply::alphaSkeleton(localOperator,ig,lfsu_s,x_s,lfsv_s,lfsu_n,x_n,lfsv_n,r_s,r_n);
  }

  // boundary integral depending on test and ansatz functions
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig,
                       const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                       R& r_s) const
  {
    Dune::PDELab::LocalOperatorApply::alphaBoundary(localOperator,ig,lfsu_s,x_s,lfsv_s,r_s);
  }

  //////////
  // lambda methods
  //////////

  // volume integral depending only on test functions
  template<typename EG, typename LFSV, typename R>
  void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
  {
    Dune::PDELab::LocalOperatorApply::lambdaVolume(localOperator,eg,lfsv,r);
  }

  template<typename EG, typename LFSV, typename R>
  void lambda_volume_post_skeleton (const EG& eg, const LFSV& lfsv, R& r) const
  {
    Dune::PDELab::LocalOperatorApply::lambdaVolumePostSkeleton(localOperator,eg,lfsv,r);
  }
  template<typename IG, typename LFSV, typename R>
  void lambda_skeleton(const IG& ig,
		       const LFSV& lfsv_s, const LFSV& lfsv_n,
		       R& r_s, R& r_n) const
  {
    Dune::PDELab::LocalOperatorApply::lambdaSkeleton(localOperator,ig,lfsv_s,lfsv_n,r_s,r_n);
  }
  template<typename IG, typename LFSV, typename R>
  void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
  {
    Dune::PDELab::LocalOperatorApply::lambdaBoundary(localOperator,ig,lfsv,r);
  }

  //////////
  // jacobian methods
  //////////
      
  template<typename EG, typename LFSU, typename X, typename LFSV, typename MAT>
  void jacobian_volume (const EG& eg,
                        const LFSU& lfsu,
                        const X& x,
                        const LFSV& lfsv,
                        MAT& mat) const
  {
    // any entity contributes to exactly one snippet!
    // volume contributes to the volume snippet this element is in
    
    // call local operator
    auto container_before = mat.container(); // copy container, so we may compute what has changed
    Dune::PDELab::LocalOperatorApply::jacobianVolume(localOperator,eg,lfsu,x,lfsv,mat);

    // get underlying container for local stiffness matrix
    const auto& container = mat.container();

    // load information about the indices
    lfs.bind(eg.entity());
    lfscache.update();
    std::vector<unsigned> localdata(lfs.size());
    view.bind(lfscache);
    view.read(localdata); // read inside element
    view.unbind();

    // scatter data
    auto volumesnippetnumber = dd.element_to_volumesnippetnumber[indexset.index(eg.entity())];
    {
      // so this element is in a skeleton snippet. Lets subtract what was computed for this volume element
      std::vector<size_t> index_in_snippet(lfs.size());
      std::vector<size_t> comp_in_snippet(lfs.size());
      for (unsigned i = 0; i<lfs.size(); ++i)
	{
	  comp_in_snippet[i] = localdata[i]%blocksize;
	  auto globali = localdata[i]/blocksize;
	  auto it = dddm.volumesnippetnumber_global_to_local[globali].find(volumesnippetnumber);
	  if (it==dddm.volumesnippetnumber_global_to_local[globali].end())
	    DUNE_THROW(Dune::Exception, "jacobian_volume: could not do global to local.");
	  index_in_snippet[i] = it->second;
	}
	
      for (size_t i = 0; i<lfs.size(); ++i)
	for (size_t j = 0; j<lfs.size(); ++j)
	  (*(*pvolume_snippet_matrices)[volumesnippetnumber])[index_in_snippet[i]][index_in_snippet[j]][comp_in_snippet[i]][comp_in_snippet[j]] +=
	      container(lfsv,i,lfsu,j) - container_before(lfsv,i,lfsu,j);

      //std::cout << "element " << indexset.index(eg.entity()) << " volume data written to volume snippet " << volumesnippetnumber << std::endl;
    }
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename MAT>
  void jacobian_volume_post_skeleton (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, MAT& mat) const
  {
    // any entity contributes to exactly one snippet!
    // volume contributes to the volume snippet this element is in

    // call local operator
    auto container_before = mat.container(); // copy container, so we may compute what has changed
    Dune::PDELab::LocalOperatorApply::jacobianVolumePostSkeleton(localOperator,eg,lfsu,x,lfsv,mat);

    // get underlying container for local stiffness matrix
    const auto& container = mat.container();

    // load information about the indices
    lfs.bind(eg.entity());
    lfscache.update();
    std::vector<unsigned> localdata(lfs.size());
    view.bind(lfscache);
    view.read(localdata); // read inside element
    view.unbind();

    // volume integral ALWAYS contributes to volume snippet and NEVER to a skeleton snippet
    // therefore we have to write it here
    auto volumesnippetnumber = dd.element_to_volumesnippetnumber[indexset.index(eg.entity())];
    {
      // so this element is in a skeleton snippet. Lets subtract what was computed for this volume element
      std::vector<size_t> index_in_snippet(lfs.size());
      std::vector<size_t> comp_in_snippet(lfs.size());
      for (unsigned i = 0; i<lfs.size(); ++i)
	{
	  comp_in_snippet[i] = localdata[i]%blocksize;
	  auto globali = localdata[i]/blocksize;
	  auto it = dddm.volumesnippetnumber_global_to_local[globali].find(volumesnippetnumber);
	  if (it==dddm.volumesnippetnumber_global_to_local[globali].end())
	    DUNE_THROW(Dune::Exception, "jacobian_volume: could not do global to local.");
	  index_in_snippet[i] = it->second;
	}
	
      for (size_t i = 0; i<lfs.size(); ++i)
	for (size_t j = 0; j<lfs.size(); ++j)
	  (*(*pvolume_snippet_matrices)[volumesnippetnumber])[index_in_snippet[i]][index_in_snippet[j]][comp_in_snippet[i]][comp_in_snippet[j]] +=
	      container(lfsv,i,lfsu,j) - container_before(lfsv,i,lfsu,j);

      //std::cout << "element volume data" << indexset.index(eg.entity()) << " written to volume snippet " << volumesnippetnumber << std::endl;
    }
  }

  template<typename IG, typename LFSU, typename X, typename LFSV, typename MAT>
  void jacobian_skeleton (const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                          MAT& mat_ss, MAT& mat_sn,
                          MAT& mat_ns, MAT& mat_nn) const
  {
    // any entity contributes to exactly one snippet!
    // for an interior face this depends on whether it is interior to a volume snippet or it is in a skeleton snippet
    
    // call local operator
    auto container_ss_before = mat_ss.container(); // copy _ss container, so we may compute what has changed, all other containers are zero
    Dune::PDELab::LocalOperatorApply::jacobianSkeleton(localOperator,ig,lfsu_s,x_s,lfsv_s,lfsu_n,x_n,lfsv_n,mat_ss,mat_sn,mat_ns,mat_nn);

    //  now we are at an interior face of the mesh. lets rad the global degrees of freedom
    // load information about the indices
    lfs.bind(ig.inside());
    lfscache.update();
    std::vector<unsigned> localdata_s(lfs.size());
    view.bind(lfscache);
    view.read(localdata_s); // read inside element
    view.unbind();

    lfs.bind(ig.outside());
    lfscache.update();
    std::vector<unsigned> localdata_n(lfs.size());
    view.bind(lfscache);
    view.read(localdata_n); // read outside element
    view.unbind();

    // get underlying container for local stiffness matrix
    const auto& container_ss = mat_ss.container();
    const auto& container_sn = mat_sn.container();
    const auto& container_ns = mat_ns.container();
    const auto& container_nn = mat_nn.container();

    // get volume snippet numbers of both elements
    auto volumesnippetnumber_s = dd.element_to_volumesnippetnumber[indexset.index(ig.inside())];
    auto volumesnippetnumber_n = dd.element_to_volumesnippetnumber[indexset.index(ig.outside())];

    // now follows the logic
    if (volumesnippetnumber_s==volumesnippetnumber_n)
      {
	// this face is interior to a volume snippet

	// prepare the global to local mapping for both elements
	std::vector<size_t> index_in_snippet_s(lfs.size());
	std::vector<size_t> comp_in_snippet_s(lfs.size());
	for (unsigned i = 0; i<lfs.size(); ++i)
	  {
	    comp_in_snippet_s[i] = localdata_s[i]%blocksize;
	    auto globali = localdata_s[i]/blocksize;
	    auto it = dddm.volumesnippetnumber_global_to_local[globali].find(volumesnippetnumber_s);
	    if (it==dddm.volumesnippetnumber_global_to_local[globali].end())
	      DUNE_THROW(Dune::Exception, "jacobian_skeleton: could not do global to local (1).");
	    index_in_snippet_s[i] = it->second;
	  }
	std::vector<size_t> index_in_snippet_n(lfs.size());
	std::vector<size_t> comp_in_snippet_n(lfs.size());
	for (unsigned i = 0; i<lfs.size(); ++i)
	  {
	    comp_in_snippet_n[i] = localdata_n[i]%blocksize;
	    auto globali = localdata_n[i]/blocksize;
	    auto it = dddm.volumesnippetnumber_global_to_local[globali].find(volumesnippetnumber_n);
	    if (it==dddm.volumesnippetnumber_global_to_local[globali].end())
	      DUNE_THROW(Dune::Exception, "jacobian_skeleton: could not do global to local (2).");
	    index_in_snippet_n[i] = it->second;
	  }

	// we write contribution to all four submatrices
	// note: volumesnippetnumber_s == volumesnippetnumber_n so we may take either
	for (size_t i = 0; i<lfs.size(); ++i)
	  for (size_t j = 0; j<lfs.size(); ++j)
	    {
	      (*(*pvolume_snippet_matrices)[volumesnippetnumber_s])[index_in_snippet_s[i]][index_in_snippet_s[j]][comp_in_snippet_s[i]][comp_in_snippet_s[j]]
		+= container_ss(lfsu_s,i,lfsu_n,j)-container_ss_before(lfsu_s,i,lfsu_s,j); 
	      (*(*pvolume_snippet_matrices)[volumesnippetnumber_s])[index_in_snippet_s[i]][index_in_snippet_n[j]][comp_in_snippet_s[i]][comp_in_snippet_n[j]] += container_sn(lfsu_s,i,lfsu_n,j); 
	      (*(*pvolume_snippet_matrices)[volumesnippetnumber_s])[index_in_snippet_n[i]][index_in_snippet_s[j]][comp_in_snippet_n[i]][comp_in_snippet_s[j]] += container_ns(lfsu_n,i,lfsu_s,j); 
	      (*(*pvolume_snippet_matrices)[volumesnippetnumber_s])[index_in_snippet_n[i]][index_in_snippet_n[j]][comp_in_snippet_n[i]][comp_in_snippet_n[j]] += container_nn(lfsu_n,i,lfsu_n,j);
	    }
	
	//std::cout << "face " << indexset.index(ig.inside()) << "|" << indexset.index(ig.outside()) << " scattered to volume snippet " << volumesnippetnumber_s << std::endl;
      }
    else
      {
	// this face is at the interface between two volume snippets
	// so it will contribute to exactly one skeleton snippet
	std::set<size_t> key; // a key to index the interface between two snippets
	key.insert(volumesnippetnumber_s);
	key.insert(volumesnippetnumber_n);
	auto skeletonsnippetnumber = dd.skeletonsnippet_to_number[key];

	// prepare the global to local mapping for both elements
	std::vector<size_t> index_in_snippet_s(lfs.size());
	std::vector<size_t> comp_in_snippet_s(lfs.size());
	for (unsigned i = 0; i<lfs.size(); ++i)
	  {
	    comp_in_snippet_s[i] = localdata_s[i]%blocksize;
	    auto globali = localdata_s[i]/blocksize;
	    auto it = dddm.skeletonsnippetnumber_global_to_local[globali].find(skeletonsnippetnumber);
	    if (it==dddm.skeletonsnippetnumber_global_to_local[globali].end())
	      DUNE_THROW(Dune::Exception, "jacobian_skeleton: could not do global to local (3).");
	    index_in_snippet_s[i] = it->second;
	  }
	std::vector<size_t> index_in_snippet_n(lfs.size());
	std::vector<size_t> comp_in_snippet_n(lfs.size());
	for (unsigned i = 0; i<lfs.size(); ++i)
	  {
	    comp_in_snippet_n[i] = localdata_n[i]%blocksize;
	    auto globali = localdata_n[i]/blocksize;
	    auto it = dddm.skeletonsnippetnumber_global_to_local[globali].find(skeletonsnippetnumber);
	    if (it==dddm.skeletonsnippetnumber_global_to_local[globali].end())
	      DUNE_THROW(Dune::Exception, "jacobian_skeleton: could not do global to local (4).");
	    index_in_snippet_n[i] = it->second;
	  }

	// we write contribution to all four submatrices
	if (pskeleton_snippet_matrices==nullptr)
	  DUNE_THROW(Dune::Exception, "jacobian_skeleton: pskeleton_snippet_matrices is null and should not be.");
	for (size_t i = 0; i<lfs.size(); ++i)
	  for (size_t j = 0; j<lfs.size(); ++j)
	    {
	      (*(*pskeleton_snippet_matrices)[skeletonsnippetnumber])[index_in_snippet_s[i]][index_in_snippet_s[j]][comp_in_snippet_s[i]][comp_in_snippet_s[j]]
		+= container_ss(lfsu_s,i,lfsu_n,j)-container_ss_before(lfsu_s,i,lfsu_s,j); 
	      (*(*pskeleton_snippet_matrices)[skeletonsnippetnumber])[index_in_snippet_s[i]][index_in_snippet_n[j]][comp_in_snippet_s[i]][comp_in_snippet_n[j]] += container_sn(lfsu_s,i,lfsu_n,j); 
	      (*(*pskeleton_snippet_matrices)[skeletonsnippetnumber])[index_in_snippet_n[i]][index_in_snippet_s[j]][comp_in_snippet_n[i]][comp_in_snippet_s[j]] += container_ns(lfsu_n,i,lfsu_s,j); 
	      (*(*pskeleton_snippet_matrices)[skeletonsnippetnumber])[index_in_snippet_n[i]][index_in_snippet_n[j]][comp_in_snippet_n[i]][comp_in_snippet_n[j]] += container_nn(lfsu_n,i,lfsu_n,j);
	    }
	
	//std::cout << "face " << indexset.index(ig.inside()) << "|" << indexset.index(ig.outside()) << " scattered to skeleton snippet " << skeletonsnippetnumber << std::endl;
      }
  }

  template<typename IG, typename LFSU, typename X, typename LFSV, typename MAT>
  void jacobian_boundary (const IG& ig,
                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                          MAT& mat_ss) const
  {
    // call local operator
    auto container_before = mat_ss.container(); // copy container, so we may compute what has changed
    Dune::PDELab::LocalOperatorApply::jacobianBoundary(localOperator,ig,lfsu_s,x_s,lfsv_s,mat_ss);

    // get underlying container for local stiffness matrix
    const auto& container = mat_ss.container();

    // load information about the indices
    lfs.bind(ig.inside());
    lfscache.update();
    std::vector<unsigned> localdata(lfs.size());
    view.bind(lfscache);
    view.read(localdata); // read inside element
    view.unbind();

    // scatter data
    auto volumesnippetnumber = dd.element_to_volumesnippetnumber[indexset.index(ig.inside())];
    {
      // so this element is in a skeleton snippet. Lets subtract what was computed for this volume element
      std::vector<size_t> index_in_snippet(lfs.size());
      std::vector<size_t> comp_in_snippet(lfs.size());
      for (unsigned i = 0; i<lfs.size(); ++i)
	{
	  comp_in_snippet[i] = localdata[i]%blocksize;
	  auto globali = localdata[i]/blocksize;
	  auto it = dddm.volumesnippetnumber_global_to_local[globali].find(volumesnippetnumber);
	  if (it==dddm.volumesnippetnumber_global_to_local[globali].end())
	    DUNE_THROW(Dune::Exception, "jacobian_volume: could not do global to local.");
	  index_in_snippet[i] = it->second;
	}
	
      for (size_t i = 0; i<lfs.size(); ++i)
	for (size_t j = 0; j<lfs.size(); ++j)
	  (*(*pvolume_snippet_matrices)[volumesnippetnumber])[index_in_snippet[i]][index_in_snippet[j]][comp_in_snippet[i]][comp_in_snippet[j]] +=
	    container(lfsv_s,i,lfsu_s,j) - container_before(lfsv_s,i,lfsu_s,j);

      //std::cout << "element " << indexset.index(ig.inside()) << " boundary data written to volume snippet " << volumesnippetnumber << std::endl;
    }
  }

  //////////
  // jacobian apply methods
  //////////

};

#endif
