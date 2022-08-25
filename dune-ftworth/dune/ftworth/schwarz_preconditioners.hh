#ifndef Udune_ftworth_schwarz_preconditioners_HH
#define Udune_ftworth_schwarz_preconditioners_HH

#include <dune/istl/umfpack.hh>

#include"matrixutilities.hh"
#include "arpack_geneo_wrapper.hh"

template<typename field_type>
std::vector<std::vector<field_type>>
get_boundary_mask (std::shared_ptr<std::vector<std::vector<size_t>>> local_to_global,
		   std::shared_ptr<std::vector<std::vector<bool>>> boundary_dofs)
{
  auto subdomains = local_to_global->size();
  std::vector<std::vector<field_type>> boundary_mask(subdomains);
  for (unsigned p=0; p<subdomains; p++)
    {
      auto n = (*local_to_global)[p].size(); // dof blocks in subdomain p
      boundary_mask[p].resize(n);
      for (unsigned i=0; i<n; i++)
	boundary_mask[p][i] = ((*boundary_dofs)[p][i]) ? 0.0 : 1.0;
    }
  return boundary_mask;
}

template<typename field_type>
std::vector<std::vector<field_type>>
get_partition_of_unity_standard (std::shared_ptr<std::vector<std::vector<size_t>>> local_to_global,
				 std::shared_ptr<std::vector<std::map<size_t,size_t>>> global_to_local,
				 std::shared_ptr<std::vector<std::vector<bool>>> boundary_dofs)
{
  // number of subdomains we have
  auto subdomains = local_to_global->size();

  // make a global vector holding the number of subdomains where dof is NOT on the boundary
  std::vector<int> k(global_to_local->size(),0);
  for (unsigned p=0; p<subdomains; p++)
    for (unsigned i=0; i<(*local_to_global)[p].size(); i++)
      if (!(*boundary_dofs)[p][i])
	k[(*local_to_global)[p][i]] += 1;

  // now compute partition of unity
  std::vector<std::vector<field_type>> pu(subdomains);
  for (unsigned p=0; p<subdomains; p++)
    {
      pu[p].resize((*local_to_global)[p].size());
      for (unsigned i=0; i<(*local_to_global)[p].size(); i++)
	pu[p][i] = ((*boundary_dofs)[p][i]) ? 0.0 : (((field_type)1.0)/k[(*local_to_global)[p][i]]);
    }

  // return
  return pu;
}

template<typename MAT>
std::vector<std::vector<typename MAT::field_type>>
get_partition_of_unity_distance (std::shared_ptr<std::vector<std::shared_ptr<MAT>>> matrices,
				 std::shared_ptr<std::vector<std::vector<size_t>>> local_to_global,
				 std::shared_ptr<std::vector<std::map<size_t,size_t>>> global_to_local,
				 std::shared_ptr<std::vector<std::vector<bool>>> boundary_dofs)
{
  // types
  using field_type = typename MAT::field_type;
  
  // number of subdomains we have
  auto subdomains = local_to_global->size();

  // compute local distance from boundary in each subdomain
  const int maxdistance = 1<<24;
  std::vector<std::vector<int>>localdistance(subdomains);
  for (unsigned p=0; p<subdomains; p++)
    {
      // dof blocks in subdomain p
      auto n = (*local_to_global)[p].size(); 

      // allocate subdomain vector and initialize with a distance never reached
      localdistance[p].resize(n,maxdistance);

      // initialize boundary degrees of freedom and interior degrees of freedom
      for (unsigned i=0; i<n; i++)
	{
	  if ((*boundary_dofs)[p][i]) localdistance[p][i] = 0; // boundary case
	  if ( (*global_to_local)[(*local_to_global)[p][i]].size()==1 ) localdistance[p][i] = -1; // interior dof case
	}
      // non boundary and non interior dofs are now at maxdistance

      // get the subdomain matrix
      MAT& A = *(*matrices)[p];

      // loop until nothing changes
      bool changed = true;
      while (changed)
	{
	  // nothing changed
	  changed = false;

	  // update for local distance
	  std::vector<int> update(localdistance[p]);

	  // compute update
	  for (unsigned i=0; i<n; i++)
	    {
	      // this dof is interior, nothing changes
	      if (localdistance[p][i]<0) continue;

	      // else we iterate over the row and compute minimum distance of neighbors and myself
	      auto cIt = A[i].begin();
	      auto cEndIt = A[i].end();
	      int dist = maxdistance;
	      for (; cIt!=cEndIt; ++cIt)
		if (localdistance[p][cIt.index()]>=0) // it is a neighbor (or myself) not in the interior
		  dist = std::min(dist,localdistance[p][cIt.index()]);

	      if (dist+1<localdistance[p][i])
		{
		  update[i] = dist+1;
		  changed = true;
		}
	    }

	  if (changed) localdistance[p] = update;
	}
      // now we have computed the distance function for the dofs in the overlap
    }

  // now compute global distance
  std::vector<int> globaldistance(global_to_local->size(),0);
  for (unsigned p=0; p<subdomains; p++)
    for (unsigned i=0; i<(*local_to_global)[p].size(); i++)
      if (localdistance[p][i]>=0)
	globaldistance[(*local_to_global)[p][i]] += localdistance[p][i];

  // now compute partition of unity
  std::vector<std::vector<field_type>> pu(subdomains);
  for (unsigned p=0; p<subdomains; p++)
    {
      pu[p].resize((*local_to_global)[p].size());
      for (unsigned i=0; i<(*local_to_global)[p].size(); i++)
	if (localdistance[p][i]<0)
	  pu[p][i] = 1.0;
	else
	  pu[p][i] = ((field_type)localdistance[p][i])/globaldistance[(*local_to_global)[p][i]];
    }

  // return
  return pu;
}

#if HAVE_SUITESPARSE_UMFPACK || DOXYGEN

/** \brief Standard single level additive Schwarz preconditioner
 *
 * \author Peter Bastian
 *
 * \tparam MAT "BCRSMatrix<FieldMatrix>" type
 * \tparam VEC "BlockVector<FieldVector>" type
 *
 */
template<typename MAT, typename VEC>
class SimpleAdditiveSchwarz : public Dune::Preconditioner<VEC,VEC>
{
  using VectorBlockType = typename VEC::block_type;
  
public:
  //! \brief The domain type of the preconditioner.
  using domain_type = VEC;
  //! \brief The range type of the preconditioner.
  using range_type = VEC;
  //! \brief The field type of the preconditioner.
  using field_type = typename VEC::field_type;

  // define the category
  virtual Dune::SolverCategory::Category category() const
  {
    return Dune::SolverCategory::sequential;
  }

  //! \brief Constructor.
  SimpleAdditiveSchwarz (std::shared_ptr<std::vector<std::shared_ptr<MAT>>> matrices_,
			 std::shared_ptr<std::vector<std::vector<size_t>>> local_to_global_,
			 std::shared_ptr<std::vector<std::vector<bool>>> boundary_dofs_
			 )
    : blocksize(VectorBlockType::dimension), matrices(matrices_), local_to_global(local_to_global_), boundary_dofs(boundary_dofs_),
      subdomains(local_to_global->size()), boundary_mask(0), subdomainsolvers(local_to_global->size())
  {
    // mask vector to null out corrections on subdomain boundary
    boundary_mask = get_boundary_mask<field_type>(local_to_global,boundary_dofs);
    
    // set up UMFPACK subdomain solvers
    for (unsigned p=0; p<subdomains; p++)
      {
	// copy subdomain matrix and assemble dirichlet conditions on subdomain boundary
	MAT& Asubdomain = *(*matrices)[p]; // modify the matrix in place
	for (size_t i=0; i<Asubdomain.N(); i++)
	  if ((*boundary_dofs)[p][i])
	    {
	      // all dofs in this block are on the boundary
	      auto cIt = Asubdomain[i].begin();
	      auto cEndIt = Asubdomain[i].end();
	      for (; cIt!=cEndIt; ++cIt)
		{
		  (*cIt) = 0.0;
		  if (cIt.index()==i)
		    for (size_t comp=0; comp<blocksize; comp++) (*cIt)[comp][comp] = 1.0;	      
		}
	    }

	// store pointer to UMFPACK solver object holding the decomposed matrix
	// hopefully it does not need the original matrix anymore ...
	subdomainsolvers[p] = std::shared_ptr<Dune::UMFPack<MAT>>(new Dune::UMFPack<MAT>(Asubdomain,false));

	(*matrices)[p] = std::shared_ptr<MAT>(nullptr); // delete subdomain matrix
      }
  }

  /*!
    \brief Prepare the preconditioner.
  */
  virtual void pre (VEC& x, VEC& b)
  {
  }

  /*!
    \brief Apply the precondioner.
  */
  virtual void apply (VEC& v, const VEC& d)
  {
    // do not need to zero out v as it is already zero
    // solve in each subdomain
    for (unsigned p=0; p<subdomains; p++)
      {
	// dof blocks in subdomain p
	auto n = (*local_to_global)[p].size();

	// set up right hand side in subdomain
	VEC dlocal(n);
	for (size_t i=0; i<n; i++)
	  {
	    dlocal[i] = d[(*local_to_global)[p][i]];
	    dlocal[i] *= boundary_mask[p][i];
	  }

	// solve subdomain problem
	Dune::InverseOperatorResult stat;
	VEC vlocal(n);
	vlocal = 0.0;
	subdomainsolvers[p]->apply(vlocal,dlocal,stat);

	// accumulate correction from subdomain
	for (size_t i=0; i<n; i++)
	  v[(*local_to_global)[p][i]] += vlocal[i];
      }
  }

  /*!
    \brief Clean up.
  */
  virtual void post (VEC& x)
  {
  }

private:
  const unsigned blocksize;
  std::shared_ptr<std::vector<std::shared_ptr<MAT>>> matrices; // the subdomain matrices
  std::shared_ptr<std::vector<std::vector<size_t>>> local_to_global; // map local index to global index in each subdomain
  std::shared_ptr<std::vector<std::vector<bool>>> boundary_dofs; // flag for each degree of freedom whether it is on the subdomain boundary
  unsigned subdomains; // the number of subdomains
  std::vector<std::vector<field_type>> boundary_mask;
  std::vector<std::shared_ptr<Dune::UMFPack<MAT>>> subdomainsolvers;
};


/** \brief Single level restricted additive Schwarz preconditioner
 *
 * \author Peter Bastian
 *
 * \tparam MAT "BCRSMatrix<FieldMatrix>" type
 * \tparam VEC "BlockVector<FieldVector>" type
 *
 */
template<typename MAT, typename VEC>
class RestrictedAdditiveSchwarz : public Dune::Preconditioner<VEC,VEC>
{
  using VectorBlockType = typename VEC::block_type;
  
public:
  //! \brief The domain type of the preconditioner.
  using domain_type = VEC;
  //! \brief The range type of the preconditioner.
  using range_type = VEC;
  //! \brief The field type of the preconditioner.
  using field_type = typename VEC::field_type;

  // define the category
  virtual Dune::SolverCategory::Category category() const
  {
    return Dune::SolverCategory::sequential;
  }

  //! \brief Constructor.
  RestrictedAdditiveSchwarz (std::shared_ptr<std::vector<std::shared_ptr<MAT>>> matrices_,
			     std::shared_ptr<std::vector<std::vector<size_t>>> local_to_global_,
			     std::shared_ptr<std::vector<std::map<size_t,size_t>>> global_to_local_,
			     std::shared_ptr<std::vector<std::vector<bool>>> boundary_dofs_,
			     std::string pum="standard")
    : blocksize(VectorBlockType::dimension), matrices(matrices_), local_to_global(local_to_global_), global_to_local(global_to_local_),
      boundary_dofs(boundary_dofs_),
      subdomains(local_to_global->size()), boundary_mask(0), partition_of_unity(0), subdomainsolvers(local_to_global->size())
  {
    // mask vector to null out corrections on subdomain boundary
    boundary_mask = get_boundary_mask<field_type>(local_to_global,boundary_dofs);
    if (pum=="distance")
      partition_of_unity = get_partition_of_unity_distance(matrices,local_to_global,global_to_local,boundary_dofs);
    else
      partition_of_unity = get_partition_of_unity_standard<field_type>(local_to_global,global_to_local,boundary_dofs);
      
    // set up UMFPACK subdomain solvers
    for (unsigned p=0; p<subdomains; p++)
      {
	// copy subdomain matrix and assemble dirichlet conditions on subdomain boundary
	MAT& Asubdomain = *(*matrices)[p]; // modify the matrix in place
	for (size_t i=0; i<Asubdomain.N(); i++)
	  if ((*boundary_dofs)[p][i])
	    {
	      // all dofs in this block are on the boundary
	      auto cIt = Asubdomain[i].begin();
	      auto cEndIt = Asubdomain[i].end();
	      for (; cIt!=cEndIt; ++cIt)
		{
		  (*cIt) = 0.0;
		  if (cIt.index()==i)
		    for (size_t comp=0; comp<blocksize; comp++) (*cIt)[comp][comp] = 1.0;	      
		}
	    }

	// store pointer to UMFPACK solver object holding the decomposed matrix
	// hopefully it does not need the original matrix anymore ...
	subdomainsolvers[p] = std::shared_ptr<Dune::UMFPack<MAT>>(new Dune::UMFPack<MAT>(Asubdomain,false));

	(*matrices)[p] = std::shared_ptr<MAT>(nullptr); // delete subdomain matrix
      }
  }

  /*!
    \brief Prepare the preconditioner.
  */
  virtual void pre (VEC& x, VEC& b)
  {
  }

  /*!
    \brief Apply the precondioner.
  */
  virtual void apply (VEC& v, const VEC& d)
  {
    // do not need to zero out v as it is already zero
    // solve in each subdomain
    for (unsigned p=0; p<subdomains; p++)
      {
	// dof blocks in subdomain p
	auto n = (*local_to_global)[p].size();

	// set up right hand side in subdomain
	VEC dlocal(n);
	for (size_t i=0; i<n; i++)
	  {
	    dlocal[i] = d[(*local_to_global)[p][i]];
	    dlocal[i] *= boundary_mask[p][i];
	  }

	// solve subdomain problem
	Dune::InverseOperatorResult stat;
	VEC vlocal(n);
	vlocal = 0.0;
	subdomainsolvers[p]->apply(vlocal,dlocal,stat);

	// accumulate correction from subdomain
	for (size_t i=0; i<n; i++)
	  v[(*local_to_global)[p][i]] += vlocal[i]*partition_of_unity[p][i];
      }
  }

  /*!
    \brief Clean up.
  */
  virtual void post (VEC& x)
  {
  }

private:
  const unsigned blocksize;
  std::shared_ptr<std::vector<std::shared_ptr<MAT>>> matrices; // the subdomain matrices
  std::shared_ptr<std::vector<std::vector<size_t>>> local_to_global; // map local index to global index in each subdomain
  std::shared_ptr<std::vector<std::map<size_t,size_t>>> global_to_local; // necessary for partition of unity
  std::shared_ptr<std::vector<std::vector<bool>>> boundary_dofs; // flag for each degree of freedom whether it is on the subdomain boundary
  unsigned subdomains; // the number of subdomains
  std::vector<std::vector<field_type>> boundary_mask;
  std::vector<std::vector<field_type>> partition_of_unity;
  std::vector<std::shared_ptr<Dune::UMFPack<MAT>>> subdomainsolvers;
};


/** \brief Two-level additive Schwarz preconditioner based on truncated input vectors
 *
 * \author Peter Bastian
 *
 * \tparam MAT a "Dune::BCRSMatrix<Dune::FieldMatrix>" type
 * \tparam VEC a "Dune::BlockVector<Dune::FieldVector>" type
 *
 */
template<typename MAT, typename VEC>
class TwoLevelAdditiveSchwarz : public Dune::Preconditioner<VEC,VEC>
{
  using VectorBlockType = typename VEC::block_type;
  using CBLOCKMAT = Dune::FieldMatrix<typename VEC::field_type,1,1>;
  using CBLOCKVEC = Dune::FieldVector<typename VEC::field_type,1>;
  using CMAT = Dune::BCRSMatrix<CBLOCKMAT>; // type for coarse level matrix
  using CVEC = Dune::BlockVector<CBLOCKVEC>; // type for coarse level vector
  
public:
  //! \brief The domain type of the preconditioner.
  using domain_type = VEC;
  //! \brief The range type of the preconditioner.
  using range_type = VEC;
  //! \brief The field type of the preconditioner.
  using field_type = typename VEC::field_type;

  // define the category
  virtual Dune::SolverCategory::Category category() const
  {
    return Dune::SolverCategory::sequential;
  }

private: // the local data of this class
  // input arguments or derived from that
  unsigned subdomains; // the number of subdomains
  const unsigned blocksize; // size of the blocks in the input vectors
  std::shared_ptr<MAT> pAglobal;
  std::shared_ptr<VEC> global_dirichlet_dofs; // (*global_dirichlet_dofs)[i][j] is >0 if component j in block i is at the global dirichlet boundary
  std::shared_ptr<std::vector<bool>> floating; // (*floating)[p] is true if subdomain p is floating
  std::shared_ptr<std::vector<std::shared_ptr<MAT>>> matrices; // the subdomain matrices
  std::shared_ptr<std::vector<std::vector<size_t>>> local_to_global; // map local index to global index in each subdomain
  std::shared_ptr<std::vector<std::map<size_t,size_t>>> global_to_local; // necessary for partition of unity
  std::shared_ptr<std::vector<std::vector<bool>>> subdomain_boundary_blocks; // flag for each degree of freedom BLOCK whether it is on the subdomain boundary
  bool twolevel; // apply coarse level correction if true (the default)

  // things computed by the constructor for later use
  unsigned n_coarse; // total number of degrees of freedom in the coarse problem
  std::vector<std::vector<field_type>> boundary_mask; // boundary_mask[p][i] is 0 if dof block i in subdomain p is at the subdomain boundary, otherwise it is 1
  std::vector<std::vector<field_type>> partition_of_unity; // partition_of_unity[p][i] is one value for dof block i in subdomain p
  std::vector<std::vector<size_t>> coarse_vector_local_to_global; // coarse_vector_local_to_global[p][i] gives the global dof number in the coarse problem of vector i in subdomain p
  std::vector<std::shared_ptr<Dune::UMFPack<MAT>>> subdomainsolvers; // the LU decomposition of the subdomain problems
  std::vector<std::vector<VEC>> coarse_basis_vectors; //  coarse_basis_vectors[p][k] gives the kth coarse grid vector in subdomain p
  std::shared_ptr<CMAT> pcmat; // the coarse level matrix
  std::shared_ptr<Dune::UMFPack<CMAT>> coarselevelsolver; // LU decomposition of coarse level matrix

public:

  //! \brief get number of basis vectors in subdomain
  size_t number_of_basis_vectors (int subdomain)
  {
    if (subdomain<0 && subdomain>=subdomains)
      DUNE_THROW(Dune::Exception, "number_of_basis_vectors: subdomain out of range.");
    return coarse_basis_vectors[subdomain].size();
  }

  //! get basis vector k in subdomain p
  VEC prolongated_basis_vector (int p, int k)
  {
    if (p<0 && p>=subdomains)
      DUNE_THROW(Dune::Exception, "number_of_basis_vectors: subdomain out of range.");

    VEC v(global_to_local->size()); // the global vector to be returned
    v = 0.0;
    auto n = (*local_to_global)[p].size();
    for (size_t i=0; i<n; i++)
      v[(*local_to_global)[p][i]] = coarse_basis_vectors[p][k][i];
    return v;
  }
  
  //! \brief Constructor. It does all the work of setting up the components
  /**
   * \param pAglobal_ shared_ptr to the global system matrix
   * \param input_coarse_vectors input_coarse_vectors_[k] points to a vector of coefficient vectors representing functions to be included in the coarse grid, e.g. constants and linears
   * \param global_dirichlet_dofs_ global_dirichlet_dofs_[i][j] is 0 if global dof block i, component j is at Dirichlet boundary, 1 else
   * \param floating_ floating_[p] is true if subdomain p touches the globakl Dirichlet boundary
   * \param matrices_ (*matrices_)[p] points to Neumann matrix in subdomain p
   * \param local_to_global_ (*local_to_global_)[p][i] maps dof block number i in subdomain p to its global number, component number is not changed and therefore not included
   * \param global_to_local_ (*global_to_local_)[i][p] maps global dof block number i to its local dof block number in subdomain p. Caution: may only be called if global dof i exists in subdomain p, otherwise use find!
   * \param subdomain_boundary_blocks_ (*subdomain_boundary_blocks_)[p][i] is true if dof block i in subdomain p is at the Dirichlet boundary
   * \param pum a string that determines which partition of unity to use. Use "standard" for 1/k_0 or "distance" for the distance based pu
   * \param twolevel_ apply coarse level correction when true
   */
  TwoLevelAdditiveSchwarz (std::shared_ptr<MAT> pAglobal_, 
			   std::vector<std::shared_ptr<VEC>> input_coarse_vectors,
			   std::shared_ptr<VEC> global_dirichlet_dofs_,
			   std::shared_ptr<std::vector<bool>> floating_,
			   std::shared_ptr<std::vector<std::shared_ptr<MAT>>> matrices_,
			   std::shared_ptr<std::vector<std::vector<size_t>>> local_to_global_,
			   std::shared_ptr<std::vector<std::map<size_t,size_t>>> global_to_local_,
			   std::shared_ptr<std::vector<std::vector<bool>>> subdomain_boundary_blocks_,
			   std::string pum,
			   bool twolevel_)
  : subdomains(local_to_global_->size()),
    blocksize(VectorBlockType::dimension),
    pAglobal(pAglobal_),
    global_dirichlet_dofs(global_dirichlet_dofs_),
    floating(floating_),
    matrices(matrices_),
    local_to_global(local_to_global_),
    global_to_local(global_to_local_),
    subdomain_boundary_blocks(subdomain_boundary_blocks_),
    twolevel(twolevel_),
    // 
    coarse_vector_local_to_global(local_to_global_->size()),
    subdomainsolvers(local_to_global_->size()),
    coarse_basis_vectors(local_to_global_->size())
  {
    std::cout << "preparing two level preconditioner ... " << std::endl;
    // compute mask vector to null out corrections on subdomain boundary
    boundary_mask = get_boundary_mask<field_type>(local_to_global,subdomain_boundary_blocks);

    // make a partition of unity
    if (pum=="distance")
      partition_of_unity = get_partition_of_unity_distance(matrices,local_to_global,global_to_local,subdomain_boundary_blocks);
    else
      partition_of_unity = get_partition_of_unity_standard<field_type>(local_to_global,global_to_local,subdomain_boundary_blocks);

    // make the coarse problem; this requires several steps

    // 1) build subdomain graph
    std::vector<std::set<size_t>> subdomain_graph(subdomains); // subdomain_graph[p] gives the subdomains having joint degrees of freedom with subdomain p, including p
    for (size_t j=0; j<global_to_local->size(); j++) // loop over global dof numbers
      for (auto itp = (*global_to_local)[j].begin(); itp!=(*global_to_local)[j].end(); ++itp) // loop of subdomains sharing this degree of freedom
	for (auto itq = (*global_to_local)[j].begin(); itq!=(*global_to_local)[j].end(); ++itq) // loop of subdomains sharing this degree of freedom
	  subdomain_graph[itp->first].insert(itq->first);
    std::cout << "subdomain graph done" << std::endl;

    // 2) create the coarse vectors in each subdomain by
    //     a) picking out entries from the global input vectors
    //     b) multiplying with the partition of unity
    //     c) nulling out dofs on the global dirchlet boundary
    // this way, a prolongation to the global Dirichlet boundary will produce a zero contribution
    //
    // this is the only step that needs to be replaced by the eigenvector computation for GenEO!
    //
    auto n_coarse_basis_vectors_per_subdomain = input_coarse_vectors.size(); // here the number of coarse vectors in each subdomain is the same
    for (unsigned p=0; p<subdomains; p++)
      {
	// create the subdomain vectors
	auto n = (*local_to_global)[p].size(); // size of the coarse vectors in subdomain p (in blocks)
	VEC vec(n); // make a template
	vec = 0.0;
	coarse_basis_vectors[p].resize(n_coarse_basis_vectors_per_subdomain,vec); // create the subdomain istl vectors; now coarse_basis_vectors[p][k] is an ISTL vector on subdomain p

	// fill them
	for (size_t k=0; k<n_coarse_basis_vectors_per_subdomain; k++)
	  {
	    // pick out values from input vector
	    for (size_t i=0; i<n; i++)
	      coarse_basis_vectors[p][k][i] = (*input_coarse_vectors[k])[(*local_to_global)[p][i]]; // pick out

	    // normalize
	    if (k==0)
	      {
		// the first is assumed to be the constant. Just rescale to norm 1
		coarse_basis_vectors[p][k] *= 1.0/coarse_basis_vectors[p][k].two_norm();
	      }
	    else
	      {
		// if it is not the first then we want to orthogonalize wrt. the constant. Let v_k be the current vector and v_0 the first (representing the constant).
		// Then seek z_k such that 
		// so 1) z_k = beta*(v_k - alpha*v_0)
		//    2) z_k*v_0 = 0 => z_k*v_0/beta = 0 = (v_k - alpha*v_0)*v_0 => v_k*v_0 - alpha*v_0*v_0 = 0 => alpha = v_k*v_0/(v_0*v_0) 
		//    3) z_k*z_k = 1
		auto alpha = coarse_basis_vectors[p][k] * coarse_basis_vectors[p][0]; // we assume that v_0 is already normalized
		coarse_basis_vectors[p][k].axpy(-alpha,coarse_basis_vectors[p][0]); // orthogonalize
		coarse_basis_vectors[p][k] *= 1.0/coarse_basis_vectors[p][k].two_norm(); // normalize
	      }
	  }

	// multiply with partition of unity and zero on Dirichlet boundary
	for (size_t k=0; k<n_coarse_basis_vectors_per_subdomain; k++)
	  for (size_t i=0; i<n; i++)
	    {
	      auto iglobal = (*local_to_global)[p][i];
	      coarse_basis_vectors[p][k][i] *= partition_of_unity[p][i]; // multiply with partition of unity, puts zero on subdomain boundary
	      for (size_t j=0; j<blocksize; j++)
		coarse_basis_vectors[p][k][i][j] *= (*global_dirichlet_dofs)[iglobal][j]; // zero out dof on global Dirichlet boundary	      
	    }
      }

    // 3) create subdomain to global mapping for the coarse problem, i.e. from subdomain and local coarse vector number to global number of coarse vectors
    // from here on this should be generic in the number of coarse vectors per subdomain
    std::vector<size_t> coarse_vector_global_to_subdomain;
    std::vector<size_t> coarse_vector_global_to_local;
    for (unsigned p=0; p<subdomains; p++)
      {
	coarse_vector_local_to_global[p].resize(coarse_basis_vectors[p].size()); // resize to number of basis vectors in this subdomain
	for (size_t k=0; k<coarse_basis_vectors[p].size(); k++) // loop over all coarse vectors in this subdomain
	  {
	    coarse_vector_global_to_subdomain.push_back(p); // this coarse vector is in subdomain p; push_back assigns the global number!
	    coarse_vector_global_to_local.push_back(k); // it has the local number k; push_back assigns the same global number!
	    coarse_vector_local_to_global[p][k] = coarse_vector_global_to_subdomain.size()-1; // could also use coarse_vector_global_to_local.size()-1
	  }
      }
    std::cout << coarse_vector_global_to_subdomain.size() << " coarse grid basis vectors" << std::endl;

    // 4) create the sparse coarse level matrix
    // This matrix in general has a different type than the input matrices since the blocksize is always 1 regardless of the blocksize of the input matrices.
    n_coarse = coarse_vector_global_to_subdomain.size(); // now we know the size of the coarse problem
    unsigned nz_coarse = 0; // compute number of nonzeroes
    for (unsigned p=0; p<subdomain_graph.size(); p++)
      for (auto q : subdomain_graph[p])
	nz_coarse += coarse_basis_vectors[p].size()*coarse_basis_vectors[q].size();
    std::cout << nz_coarse << " nonzeroes in coarse level matrix" << std::endl;
    pcmat = std::shared_ptr<CMAT>(new CMAT(n_coarse,n_coarse,nz_coarse,CMAT::row_wise)); // allocate matrix in rowwise creation mode
    for (auto row=pcmat->createbegin(); row!=pcmat->createend(); ++row) // fill nonzeroes
      {
	auto i = row.index(); // this is the global row index
	auto p = coarse_vector_global_to_subdomain[i]; // the subdomain this coarse vector comes from
	for (auto q : subdomain_graph[p]) // loop over neighboring subdomains
	  for (size_t k=0; k<coarse_basis_vectors[q].size(); k++) // loop over coarse vectors in neighboring subdomain
	    {
	      auto j = coarse_vector_local_to_global[q][k]; // the global number of this coarse vector
	      row.insert(j);
	    }
      }
    std::cout << "coarse level matrix created" << std::endl;
 
    // 5) compute entries of coarse level matrix
    for (unsigned prow=0; prow<subdomains; prow++) // loop over all subdomains (row index)
      {
	// copy the subdomain matrix and fill it with entries of the global matrix
	// this is needed to get the correct entries on the boundary of the subdomain
	MAT Asubdomain(*(*matrices)[prow]);
	for (size_t i=0; i<Asubdomain.N(); i++) // loop over rows
	  {
	    auto iglobal = (*local_to_global)[prow][i];
	    auto cIt = Asubdomain[i].begin(); 
	    auto cEndIt = Asubdomain[i].end();
	    for (; cIt!=cEndIt; ++cIt) // loop over columns
	      {
		auto j = cIt.index();
		auto jglobal = (*local_to_global)[prow][j];
		(*cIt) = (*pAglobal)[iglobal][jglobal];
	      }
	  }
	std::cout << prow << ": subdomain matrix in coarse system computation is symmetric: " << is_symmetric(Asubdomain," ",true,8e-12) << std::endl;
	
	// now do the triple matrix products
	VEC v1(Asubdomain.N()); // two vectors on the subdomain to compute scalar products
	VEC v2(Asubdomain.N());
	for (size_t krow=0; krow<coarse_basis_vectors[prow].size(); krow++) // loop over all coarse grid functions in this subdomain
	  {
	    // row number in coarse system
	    auto icoarse = coarse_vector_local_to_global[prow][krow];
	    
	    // matrix vector product
	    Asubdomain.mv(coarse_basis_vectors[prow][krow],v1);

	    // now do all the scalar products with the other basis vectors
	    for (auto pcol : subdomain_graph[prow]) // loop over all neighboring subdomains (including myself)
	      for (size_t kcol=0; kcol<coarse_basis_vectors[pcol].size(); kcol++) // loop over all coarse grid functions in the neighboring subdomain
		{
		  // column number in coarse system
		  auto jcoarse = coarse_vector_local_to_global[pcol][kcol];

		  // produce restriction of basis vector from subdomain pcol to subdomain prow
		  v2 = 0.0; // zero out
		  for (size_t i=0; i<Asubdomain.N(); i++)
		    {
		      auto iglobal = (*local_to_global)[prow][i]; // global dof block for our local index
		      auto it = (*global_to_local)[iglobal].find(pcol); // test if this global dof block exists in neighboring subdomain
		      if (it!=(*global_to_local)[iglobal].end()) // if yes, we have a common entry
			{
			  auto index_in_other_subdomain = it->second;
			  v2[i] = coarse_basis_vectors[pcol][kcol][index_in_other_subdomain]; // restriction of coarse vector in subdomain pcol in subdomain prow
			}
		    }

		  // now we have a matrix entry
		  (*pcmat)[icoarse][jcoarse] = v2*v1;
		}
	  }
      }
    std::cout << "computed entries of coarse system" << std::endl;
    //Dune::printmatrix(std::cout, *pcmat, "coarse level matrix", "");

    // // write coarse matrix to file
    // {
    //   std::ostringstream s1;
    //   s1 << "DGcoarsematrix_" << pcmat->M() << "_" << subdomains << ".txt";
    //   Dune::writeMatrixToMatlab(*pcmat,s1.str());
    // }

    // check symmetry
    std::cout << "coarse level matrix is symmetric: " << is_symmetric(*pcmat," ",true,8e-12) << std::endl;

    // now we may compute the LU decomposition of the coarse system
    std::cout << "LU decomposition of coarse system" << std::endl;
    coarselevelsolver = std::shared_ptr<Dune::UMFPack<CMAT>>(new Dune::UMFPack<CMAT>(*pcmat,false));
    std::cout << "... done" << std::endl;
    
    // set up UMFPACK subdomain solvers
    std::cout << "LU decomposition of subdomain systems" << std::endl;
    for (unsigned p=0; p<subdomains; p++)
      {
	// copy subdomain matrix and assemble dirichlet conditions on subdomain boundary
	MAT& Asubdomain = *(*matrices)[p]; // modify the matrix in place as we do not need it anymore
	for (size_t i=0; i<Asubdomain.N(); i++)
	  if ((*subdomain_boundary_blocks)[p][i])
	    {
	      // all dofs in this block are on the boundary
	      auto cIt = Asubdomain[i].begin();
	      auto cEndIt = Asubdomain[i].end();
	      for (; cIt!=cEndIt; ++cIt)
		{
		  (*cIt) = 0.0;
		  if (cIt.index()==i)
		    for (size_t comp=0; comp<blocksize; comp++) (*cIt)[comp][comp] = 1.0;	      
		}
	    }

	// store pointer to UMFPACK solver object holding the decomposed matrix
	// hopefully it does not need the original matrix anymore ...
	subdomainsolvers[p] = std::shared_ptr<Dune::UMFPack<MAT>>(new Dune::UMFPack<MAT>(Asubdomain,false));

	(*matrices)[p] = std::shared_ptr<MAT>(nullptr); // delete subdomain matrix
      }
    std::cout << "... done" << std::endl;
    std::cout << "coarse level correction: " << twolevel << std::endl;
    std::cout << "preconditioner ready." << std::endl;
  }

  /*!
    \brief Prepare the preconditioner.
  */
  virtual void pre (VEC& x, VEC& b)
  {
  }

  /*!
    \brief Apply the precondioner.
  */
  virtual void apply (VEC& v, const VEC& d)
  {
    // do not need to zero out v as it is already zero
    // solve in each subdomain
    for (unsigned p=0; p<subdomains; p++)
      {
	// dof blocks in subdomain p
	auto n = (*local_to_global)[p].size();

	// set up right hand side in subdomain
	VEC dlocal(n);
	for (size_t i=0; i<n; i++)
	  {
	    dlocal[i] = d[(*local_to_global)[p][i]];
	    dlocal[i] *= boundary_mask[p][i];
	  }

	// solve subdomain problem
	Dune::InverseOperatorResult stat;
	VEC vlocal(n);
	vlocal = 0.0;
	subdomainsolvers[p]->apply(vlocal,dlocal,stat);

	// accumulate correction from subdomain
	for (size_t i=0; i<n; i++)
	  v[(*local_to_global)[p][i]] += vlocal[i]; //*partition_of_unity[p][i];
      }

    // check if coarse level correction is requested
    if (!twolevel) return;

    // prepare and solve coarse problem
    // the coarse grid defect
    CVEC dcoarse(n_coarse);
    dcoarse = 0.0;

    // restriction
    for (unsigned p=0; p<subdomains; p++)
      for (size_t k=0; k<coarse_basis_vectors[p].size(); k++)
	{
	  auto icoarse = coarse_vector_local_to_global[p][k];
	  auto n = (*local_to_global)[p].size();
	  for (size_t i=0; i<n; i++)
	    dcoarse[icoarse] += d[(*local_to_global)[p][i]]*coarse_basis_vectors[p][k][i];
	}
	  
    // solve subdomain problem
    Dune::InverseOperatorResult stat;
    CVEC vcoarse(n_coarse);
    vcoarse = 0.0;
    coarselevelsolver->apply(vcoarse,dcoarse,stat);

    // prolongation
    for (unsigned p=0; p<subdomains; p++)
      for (size_t k=0; k<coarse_basis_vectors[p].size(); k++)
	{
	  auto icoarse = coarse_vector_local_to_global[p][k];
	  auto n = (*local_to_global)[p].size();
	  for (size_t i=0; i<n; i++)
	    v[(*local_to_global)[p][i]] += coarse_basis_vectors[p][k][i]*vcoarse[icoarse][0];
	}
  }

  /*!
    \brief Clean up.
  */
  virtual void post (VEC& x)
  {
  }

};



/** \brief Two-level additive Schwarz preconditioner based on truncated input vectors
 *
 * \author Peter Bastian
 *
 * \tparam MAT a "Dune::BCRSMatrix<Dune::FieldMatrix>" type
 * \tparam VEC a "Dune::BlockVector<Dune::FieldVector>" type
 *
 */
template<typename MAT, typename VEC>
class PetersTwoLevelGenEO : public Dune::Preconditioner<VEC,VEC>
{
  using VectorBlockType = typename VEC::block_type;
  using CBLOCKMAT = Dune::FieldMatrix<typename VEC::field_type,1,1>;
  using CBLOCKVEC = Dune::FieldVector<typename VEC::field_type,1>;
  using CMAT = Dune::BCRSMatrix<CBLOCKMAT>; // type for coarse level matrix
  using CVEC = Dune::BlockVector<CBLOCKVEC>; // type for coarse level vector
  
public:
  //! \brief The domain type of the preconditioner.
  using domain_type = VEC;
  //! \brief The range type of the preconditioner.
  using range_type = VEC;
  //! \brief The field type of the preconditioner.
  using field_type = typename VEC::field_type;

  // define the category
  virtual Dune::SolverCategory::Category category() const
  {
    return Dune::SolverCategory::sequential;
  }

private: // the local data of this class
  // input arguments or derived from that
  unsigned subdomains; // the number of subdomains
  const unsigned blocksize; // size of the blocks in the input vectors
  std::shared_ptr<MAT> pAglobal;
  std::shared_ptr<VEC> global_dirichlet_dofs; // (*global_dirichlet_dofs)[i][j] is >0 if component j in block i is at the global dirichlet boundary
  std::shared_ptr<std::vector<bool>> floating; // (*floating)[p] is true if subdomain p is floating
  std::shared_ptr<std::vector<std::shared_ptr<MAT>>> matrices; // the subdomain matrices
  std::shared_ptr<std::vector<std::vector<size_t>>> local_to_global; // map local index to global index in each subdomain
  std::shared_ptr<std::vector<std::map<size_t,size_t>>> global_to_local; // necessary for partition of unity
  std::shared_ptr<std::vector<std::vector<bool>>> subdomain_boundary_blocks; // flag for each degree of freedom BLOCK whether it is on the subdomain boundary
  bool twolevel; // apply coarse level correction if true (the default)

  // things computed by the constructor for later use
  unsigned n_coarse; // total number of degrees of freedom in the coarse problem
  std::vector<std::vector<field_type>> boundary_mask; // boundary_mask[p][i] is 0 if dof block i in subdomain p is at the subdomain boundary, otherwise it is 1
  std::vector<std::vector<field_type>> partition_of_unity; // partition_of_unity[p][i] is one value for dof block i in subdomain p
  std::vector<std::vector<size_t>> coarse_vector_local_to_global; // coarse_vector_local_to_global[p][i] gives the global dof number in the coarse problem of vector i in subdomain p
  std::vector<std::shared_ptr<Dune::UMFPack<MAT>>> subdomainsolvers; // the LU decomposition of the subdomain problems
  std::vector<std::vector<VEC>> coarse_basis_vectors; //  coarse_basis_vectors[p][k] gives the kth coarse grid vector in subdomain p
  std::shared_ptr<CMAT> pcmat; // the coarse level matrix
  std::shared_ptr<Dune::UMFPack<CMAT>> coarselevelsolver; // LU decomposition of coarse level matrix

public:

  //! \brief get number of basis vectors in subdomain
  size_t number_of_basis_vectors (int subdomain)
  {
    if (subdomain<0 && subdomain>=subdomains)
      DUNE_THROW(Dune::Exception, "number_of_basis_vectors: subdomain out of range.");
    return coarse_basis_vectors[subdomain].size();
  }

  //! get basis vector k in subdomain p
  VEC prolongated_basis_vector (int p, int k)
  {
    if (p<0 && p>=subdomains)
      DUNE_THROW(Dune::Exception, "number_of_basis_vectors: subdomain out of range.");

    VEC v(global_to_local->size()); // the global vector to be returned
    v = 0.0;
    auto n = (*local_to_global)[p].size();
    for (size_t i=0; i<n; i++)
      v[(*local_to_global)[p][i]] = coarse_basis_vectors[p][k][i];
    return v;
  }
  
  //! \brief Constructor. It does all the work of setting up the components
  /**
   * \param pAglobal_ shared_ptr to the global system matrix
   * \param input_coarse_vectors input_coarse_vectors_[k] points to a vector of coefficient vectors representing functions to be included in the coarse grid, e.g. constants and linears
   * \param global_dirichlet_dofs_ global_dirichlet_dofs_[i][j] is 0 if global dof block i, component j is at Dirichlet boundary, 1 else
   * \param floating_ floating_[p] is true if subdomain p touches the globakl Dirichlet boundary
   * \param matrices_ (*matrices_)[p] points to Neumann matrix in subdomain p
   * \param local_to_global_ (*local_to_global_)[p][i] maps dof block number i in subdomain p to its global number, component number is not changed and therefore not included
   * \param global_to_local_ (*global_to_local_)[i][p] maps global dof block number i to its local dof block number in subdomain p. Caution: may only be called if global dof i exists in subdomain p, otherwise use find!
   * \param subdomain_boundary_blocks_ (*subdomain_boundary_blocks_)[p][i] is true if dof block i in subdomain p is at the Dirichlet boundary
   * \param pum a string that determines which partition of unity to use. Use "standard" for 1/k_0 or "distance" for the distance based pu
   * \param twolevel_ apply coarse level correction when true
   */
  PetersTwoLevelGenEO (std::shared_ptr<MAT> pAglobal_, 
			   std::vector<std::shared_ptr<VEC>> input_coarse_vectors,
			   std::shared_ptr<VEC> global_dirichlet_dofs_,
			   std::shared_ptr<std::vector<bool>> floating_,
			   std::shared_ptr<std::vector<std::shared_ptr<MAT>>> matrices_,
			   std::shared_ptr<std::vector<std::vector<size_t>>> local_to_global_,
			   std::shared_ptr<std::vector<std::map<size_t,size_t>>> global_to_local_,
			   std::shared_ptr<std::vector<std::vector<bool>>> subdomain_boundary_blocks_,
			   std::string pum,
			   bool twolevel_)
  : subdomains(local_to_global_->size()),
    blocksize(VectorBlockType::dimension),
    pAglobal(pAglobal_),
    global_dirichlet_dofs(global_dirichlet_dofs_),
    floating(floating_),
    matrices(matrices_),
    local_to_global(local_to_global_),
    global_to_local(global_to_local_),
    subdomain_boundary_blocks(subdomain_boundary_blocks_),
    twolevel(twolevel_),
    // 
    coarse_vector_local_to_global(local_to_global_->size()),
    subdomainsolvers(local_to_global_->size()),
    coarse_basis_vectors(local_to_global_->size())
  {
    std::cout << "preparing two level preconditioner ... " << std::endl;
    // compute mask vector to null out corrections on subdomain boundary
    boundary_mask = get_boundary_mask<field_type>(local_to_global,subdomain_boundary_blocks);

    // make a partition of unity
    if (pum=="distance")
      partition_of_unity = get_partition_of_unity_distance(matrices,local_to_global,global_to_local,subdomain_boundary_blocks);
    else
      partition_of_unity = get_partition_of_unity_standard<field_type>(local_to_global,global_to_local,subdomain_boundary_blocks);

    // make the coarse problem; this requires several steps

    // 1) build subdomain graph
    std::vector<std::set<size_t>> subdomain_graph(subdomains); // subdomain_graph[p] gives the subdomains having joint degrees of freedom with subdomain p, including p
    for (size_t j=0; j<global_to_local->size(); j++) // loop over global dof numbers
      for (auto itp = (*global_to_local)[j].begin(); itp!=(*global_to_local)[j].end(); ++itp) // loop of subdomains sharing this degree of freedom
	for (auto itq = (*global_to_local)[j].begin(); itq!=(*global_to_local)[j].end(); ++itq) // loop of subdomains sharing this degree of freedom
	  subdomain_graph[itp->first].insert(itq->first);
    std::cout << "subdomain graph done" << std::endl;

    // 2) create the coarse vectors in each subdomain by
    //     a) picking out entries from the global input vectors
    //     b) multiplying with the partition of unity
    //     c) nulling out dofs on the global dirchlet boundary
    // this way, a prolongation to the global Dirichlet boundary will produce a zero contribution
    //
    // this is the only step that needs to be replaced by the eigenvector computation for GenEO!
    //
    auto n_coarse_basis_vectors_per_subdomain = input_coarse_vectors.size(); // here the number of coarse vectors in each subdomain is the same
    for (unsigned p=0; p<subdomains; p++)
      {
	// solve the eigenvalue problem
	// A is the neumann matrix. It will be the left-hand side of our eigenproblem. (This is redundant but writing 'A' is nicer.)
	MAT& A = *(*matrices)[p]; 

	// B will be the right-hand side of our eigenproblem. It starts as a physical copy of the neumann matrix.
	MAT B(A);

	// Now multiply B from left and right with the partition of unity.
	for (auto row_iter = B.begin(); row_iter != B.end(); ++row_iter) 
	  for (auto col_iter = row_iter->begin(); col_iter != row_iter->end(); ++col_iter) 
	    *col_iter *= partition_of_unity[p][row_iter.index()] * partition_of_unity[p][col_iter.index()];

	// Setup ArpackGenEO wrapper, solving a generalized eigenproblem with lhs A.
	ArpackMLGeneo::ArPackPlusPlus_Algorithms<MAT,VEC> arpack(A);

	// set some parameters. Some or all of these should be input parameters.
	double tolerance = 0.0; // tolerance for arpack algorithm
	int number_of_eigenvalues = 20; // how many eigenvalues does arpack compute
	double shift = 0.001; // shift used in arpack algorihm
	double threshold = -1; // eigenvalues below this are considered interesting. Their eigenvectors are added to coarse_basis_vectors
	int verbose = 1;

	// make seperate vectors for arpack
	auto n = (*local_to_global)[p].size(); // size of the coarse vectors in subdomain p (in blocks)
	VEC vec(n); // make a template
	vec = 0.0;
	std::vector<VEC> eigenvectors(number_of_eigenvalues,vec); // create vector of ISTL vectors for Arpack to compute        

	// make vectors to store eigenvalues. Eigenvectors will be written directly into coarse_basis_vectors.
	std::vector<double> eigenvalues(number_of_eigenvalues,0.0);

	//====================
	bool use_geneo = true;
	//====================

	// solve GenEO eigenproblems
	if (use_geneo)
	  {
	    std::cout << "subdomain " << p << ": solve eigenvalue problem" << std::endl;
	    arpack.computeStdNonSymMinMagnitude(B,tolerance,eigenvectors,eigenvalues,shift);
	    for (size_t i=0; i<eigenvalues.size(); ++i)
	      std::cout << "  lambda_" << i<< " = " << eigenvalues[i] << std::endl;
	  }
	
	// create the subdomain vectors
	coarse_basis_vectors[p].resize(n_coarse_basis_vectors_per_subdomain,vec); // create the subdomain istl vectors; now coarse_basis_vectors[p][k] is an ISTL vector on subdomain p

	// fill them
	for (size_t k=0; k<n_coarse_basis_vectors_per_subdomain; k++)
	  {
	    // pick out values from input vector
	    if (use_geneo && number_of_eigenvalues>=n_coarse_basis_vectors_per_subdomain)
	      {
		coarse_basis_vectors[p][k] = eigenvectors[k]; // copy ev from arpack output
	      }
	    else
	      {
		for (size_t i=0; i<n; i++)
		  coarse_basis_vectors[p][k][i] = (*input_coarse_vectors[k])[(*local_to_global)[p][i]]; // pick out
	      }

	    // normalize
	    coarse_basis_vectors[p][k] *= 1.0/coarse_basis_vectors[p][k].two_norm();
	    // if (k==0)
	    //   {
	    // 	// the first is assumed to be the constant. Just rescale to norm 1
	    // 	coarse_basis_vectors[p][k] *= 1.0/coarse_basis_vectors[p][k].two_norm();
	    //   }
	    // else
	    //   {
	    // 	// if it is not the first then we want to orthogonalize wrt. the constant. Let v_k be the current vector and v_0 the first (representing the constant).
	    // 	// Then seek z_k such that 
	    // 	// so 1) z_k = beta*(v_k - alpha*v_0)
	    // 	//    2) z_k*v_0 = 0 => z_k*v_0/beta = 0 = (v_k - alpha*v_0)*v_0 => v_k*v_0 - alpha*v_0*v_0 = 0 => alpha = v_k*v_0/(v_0*v_0) 
	    // 	//    3) z_k*z_k = 1
	    // 	auto alpha = coarse_basis_vectors[p][k] * coarse_basis_vectors[p][0]; // we assume that v_0 is already normalized
	    // 	coarse_basis_vectors[p][k].axpy(-alpha,coarse_basis_vectors[p][0]); // orthogonalize
	    // 	coarse_basis_vectors[p][k] *= 1.0/coarse_basis_vectors[p][k].two_norm(); // normalize
	    //   }
	  }

	// multiply with partition of unity and zero on Dirichlet boundary
	for (size_t k=0; k<n_coarse_basis_vectors_per_subdomain; k++)
	  for (size_t i=0; i<n; i++)
	    {
	      auto iglobal = (*local_to_global)[p][i];
	      coarse_basis_vectors[p][k][i] *= partition_of_unity[p][i]; // multiply with partition of unity, puts zero on subdomain boundary
	      for (size_t j=0; j<blocksize; j++)
		coarse_basis_vectors[p][k][i][j] *= (*global_dirichlet_dofs)[iglobal][j]; // zero out dof on global Dirichlet boundary	      
	    }
      }

    // 3) create subdomain to global mapping for the coarse problem, i.e. from subdomain and local coarse vector number to global number of coarse vectors
    // from here on this should be generic in the number of coarse vectors per subdomain
    std::vector<size_t> coarse_vector_global_to_subdomain;
    std::vector<size_t> coarse_vector_global_to_local;
    for (unsigned p=0; p<subdomains; p++)
      {
	coarse_vector_local_to_global[p].resize(coarse_basis_vectors[p].size()); // resize to number of basis vectors in this subdomain
	for (size_t k=0; k<coarse_basis_vectors[p].size(); k++) // loop over all coarse vectors in this subdomain
	  {
	    coarse_vector_global_to_subdomain.push_back(p); // this coarse vector is in subdomain p; push_back assigns the global number!
	    coarse_vector_global_to_local.push_back(k); // it has the local number k; push_back assigns the same global number!
	    coarse_vector_local_to_global[p][k] = coarse_vector_global_to_subdomain.size()-1; // could also use coarse_vector_global_to_local.size()-1
	  }
      }
    std::cout << coarse_vector_global_to_subdomain.size() << " coarse grid basis vectors" << std::endl;

    // 4) create the sparse coarse level matrix
    // This matrix in general has a different type than the input matrices since the blocksize is always 1 regardless of the blocksize of the input matrices.
    n_coarse = coarse_vector_global_to_subdomain.size(); // now we know the size of the coarse problem
    unsigned nz_coarse = 0; // compute number of nonzeroes
    for (unsigned p=0; p<subdomain_graph.size(); p++)
      for (auto q : subdomain_graph[p])
	nz_coarse += coarse_basis_vectors[p].size()*coarse_basis_vectors[q].size();
    std::cout << nz_coarse << " nonzeroes in coarse level matrix" << std::endl;
    pcmat = std::shared_ptr<CMAT>(new CMAT(n_coarse,n_coarse,nz_coarse,CMAT::row_wise)); // allocate matrix in rowwise creation mode
    for (auto row=pcmat->createbegin(); row!=pcmat->createend(); ++row) // fill nonzeroes
      {
	auto i = row.index(); // this is the global row index
	auto p = coarse_vector_global_to_subdomain[i]; // the subdomain this coarse vector comes from
	for (auto q : subdomain_graph[p]) // loop over neighboring subdomains
	  for (size_t k=0; k<coarse_basis_vectors[q].size(); k++) // loop over coarse vectors in neighboring subdomain
	    {
	      auto j = coarse_vector_local_to_global[q][k]; // the global number of this coarse vector
	      row.insert(j);
	    }
      }
    std::cout << "coarse level matrix created" << std::endl;
 
    // 5) compute entries of coarse level matrix
    for (unsigned prow=0; prow<subdomains; prow++) // loop over all subdomains (row index)
      {
	// copy the subdomain matrix and fill it with entries of the global matrix
	// this is needed to get the correct entries on the boundary of the subdomain
	MAT Asubdomain(*(*matrices)[prow]);
	for (size_t i=0; i<Asubdomain.N(); i++) // loop over rows
	  {
	    auto iglobal = (*local_to_global)[prow][i];
	    auto cIt = Asubdomain[i].begin(); 
	    auto cEndIt = Asubdomain[i].end();
	    for (; cIt!=cEndIt; ++cIt) // loop over columns
	      {
		auto j = cIt.index();
		auto jglobal = (*local_to_global)[prow][j];
		(*cIt) = (*pAglobal)[iglobal][jglobal];
	      }
	  }
	std::cout << prow << ": subdomain matrix in coarse system computation is symmetric: " << is_symmetric(Asubdomain," ",true,8e-12) << std::endl;
	
	// now do the triple matrix products
	VEC v1(Asubdomain.N()); // two vectors on the subdomain to compute scalar products
	VEC v2(Asubdomain.N());
	for (size_t krow=0; krow<coarse_basis_vectors[prow].size(); krow++) // loop over all coarse grid functions in this subdomain
	  {
	    // row number in coarse system
	    auto icoarse = coarse_vector_local_to_global[prow][krow];
	    
	    // matrix vector product
	    Asubdomain.mv(coarse_basis_vectors[prow][krow],v1);

	    // now do all the scalar products with the other basis vectors
	    for (auto pcol : subdomain_graph[prow]) // loop over all neighboring subdomains (including myself)
	      for (size_t kcol=0; kcol<coarse_basis_vectors[pcol].size(); kcol++) // loop over all coarse grid functions in the neighboring subdomain
		{
		  // column number in coarse system
		  auto jcoarse = coarse_vector_local_to_global[pcol][kcol];

		  // produce restriction of basis vector from subdomain pcol to subdomain prow
		  v2 = 0.0; // zero out
		  for (size_t i=0; i<Asubdomain.N(); i++)
		    {
		      auto iglobal = (*local_to_global)[prow][i]; // global dof block for our local index
		      auto it = (*global_to_local)[iglobal].find(pcol); // test if this global dof block exists in neighboring subdomain
		      if (it!=(*global_to_local)[iglobal].end()) // if yes, we have a common entry
			{
			  auto index_in_other_subdomain = it->second;
			  v2[i] = coarse_basis_vectors[pcol][kcol][index_in_other_subdomain]; // restriction of coarse vector in subdomain pcol in subdomain prow
			}
		    }

		  // now we have a matrix entry
		  (*pcmat)[icoarse][jcoarse] = v2*v1;
		}
	  }
      }
    std::cout << "computed entries of coarse system" << std::endl;

    //Dune::printmatrix(std::cout, *pcmat, "coarse level matrix", "");

    // // write coarse matrix to file
    // {
    //   std::ostringstream s1;
    //   s1 << "DGcoarsematrix_" << pcmat->M() << "_" << subdomains << ".txt";
    //   Dune::writeMatrixToMatlab(*pcmat,s1.str());
    // }

    // check symmetry
    std::cout << "coarse level matrix is symmetric: " << is_symmetric(*pcmat," ",true,8e-12) << std::endl;

    // now we may compute the LU decomposition of the coarse system
    std::cout << "LU decomposition of coarse system" << std::endl;
    coarselevelsolver = std::shared_ptr<Dune::UMFPack<CMAT>>(new Dune::UMFPack<CMAT>(*pcmat,false));
    std::cout << "... done" << std::endl;
    
    // set up UMFPACK subdomain solvers
    std::cout << "LU decomposition of subdomain systems" << std::endl;
    for (unsigned p=0; p<subdomains; p++)
      {
	// copy subdomain matrix and assemble dirichlet conditions on subdomain boundary
	MAT& Asubdomain = *(*matrices)[p]; // modify the matrix in place as we do not need it anymore
	for (size_t i=0; i<Asubdomain.N(); i++)
	  if ((*subdomain_boundary_blocks)[p][i])
	    {
	      // all dofs in this block are on the boundary
	      auto cIt = Asubdomain[i].begin();
	      auto cEndIt = Asubdomain[i].end();
	      for (; cIt!=cEndIt; ++cIt)
		{
		  (*cIt) = 0.0;
		  if (cIt.index()==i)
		    for (size_t comp=0; comp<blocksize; comp++) (*cIt)[comp][comp] = 1.0;	      
		}
	    }

	// store pointer to UMFPACK solver object holding the decomposed matrix
	// hopefully it does not need the original matrix anymore ...
	subdomainsolvers[p] = std::shared_ptr<Dune::UMFPack<MAT>>(new Dune::UMFPack<MAT>(Asubdomain,false));

	(*matrices)[p] = std::shared_ptr<MAT>(nullptr); // delete subdomain matrix
      }
    std::cout << "... done" << std::endl;
    std::cout << "coarse level correction: " << twolevel << std::endl;
    std::cout << "preconditioner ready." << std::endl;
  }

  /*!
    \brief Prepare the preconditioner.
  */
  virtual void pre (VEC& x, VEC& b)
  {
  }

  /*!
    \brief Apply the precondioner.
  */
  virtual void apply (VEC& v, const VEC& d)
  {
    // do not need to zero out v as it is already zero
    // solve in each subdomain
    for (unsigned p=0; p<subdomains; p++)
      {
	// dof blocks in subdomain p
	auto n = (*local_to_global)[p].size();

	// set up right hand side in subdomain
	VEC dlocal(n);
	for (size_t i=0; i<n; i++)
	  {
	    dlocal[i] = d[(*local_to_global)[p][i]];
	    dlocal[i] *= boundary_mask[p][i];
	  }

	// solve subdomain problem
	Dune::InverseOperatorResult stat;
	VEC vlocal(n);
	vlocal = 0.0;
	subdomainsolvers[p]->apply(vlocal,dlocal,stat);

	// accumulate correction from subdomain
	for (size_t i=0; i<n; i++)
	  v[(*local_to_global)[p][i]] += vlocal[i]; //*partition_of_unity[p][i];
      }

    // check if coarse level correction is requested
    if (!twolevel) return;

    // prepare and solve coarse problem
    // the coarse grid defect
    CVEC dcoarse(n_coarse);
    dcoarse = 0.0;

    // restriction
    for (unsigned p=0; p<subdomains; p++)
      for (size_t k=0; k<coarse_basis_vectors[p].size(); k++)
	{
	  auto icoarse = coarse_vector_local_to_global[p][k];
	  auto n = (*local_to_global)[p].size();
	  for (size_t i=0; i<n; i++)
	    dcoarse[icoarse] += d[(*local_to_global)[p][i]]*coarse_basis_vectors[p][k][i];
	}
	  
    // solve subdomain problem
    Dune::InverseOperatorResult stat;
    CVEC vcoarse(n_coarse);
    vcoarse = 0.0;
    coarselevelsolver->apply(vcoarse,dcoarse,stat);

    // prolongation
    for (unsigned p=0; p<subdomains; p++)
      for (size_t k=0; k<coarse_basis_vectors[p].size(); k++)
	{
	  auto icoarse = coarse_vector_local_to_global[p][k];
	  auto n = (*local_to_global)[p].size();
	  for (size_t i=0; i<n; i++)
	    v[(*local_to_global)[p][i]] += coarse_basis_vectors[p][k][i]*vcoarse[icoarse][0];
	}
  }

  /*!
    \brief Clean up.
  */
  virtual void post (VEC& x)
  {
  }

};

#endif // HAVE_SUITESPARSE_UMFPACK || DOXYGEN

#endif // Udune_ftworth_HH
