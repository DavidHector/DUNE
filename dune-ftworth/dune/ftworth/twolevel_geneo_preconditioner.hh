#ifndef Udune_ftworth_twolevel_geneo_preconditioner_HH
#define Udune_ftworth_twolevel_geneo_preconditioner_HH

#include <dune/istl/umfpack.hh>
#include<dune/istl/matrixmarket.hh>

#include "schwarz_preconditioners.hh"
// #include "arpack_geneo_wrapper.hh"
#include <dune/pdelab/backend/istl/geneo/arpackpp_geneo.hh> // needs pdelab branch "geneo_nonovlp"

#if HAVE_SUITESPARSE_UMFPACK || DOXYGEN


/** \brief Two-level GenEO preconditioner
 *
 * \author Arne Strehlow
 *
 * \tparam MAT a "Dune::BCRSMatrix<Dune::FieldMatrix>" type
 * \tparam VEC a "Dune::BlockVector<Dune::FieldVector>" type
 *
 */
template<typename MAT, typename VEC>
class TwoLevelGenEO : public Dune::Preconditioner<VEC,VEC>
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
   * \param debugx vectors used for debugging
   * \param ptree Parameter tree containing input parameters from *.ini file
   */
  TwoLevelGenEO (std::shared_ptr<MAT> pAglobal_, 
		 std::vector<std::shared_ptr<VEC>> input_coarse_vectors,
		 std::shared_ptr<VEC> global_dirichlet_dofs_,
		 std::shared_ptr<std::vector<bool>> floating_,
		 std::shared_ptr<std::vector<std::shared_ptr<MAT>>> matrices_,
		 std::shared_ptr<std::vector<std::vector<size_t>>> local_to_global_,
		 std::shared_ptr<std::vector<std::map<size_t,size_t>>> global_to_local_,
		 std::shared_ptr<std::vector<std::vector<bool>>> subdomain_boundary_blocks_,
		 std::string pum,
		 bool twolevel_,
		 VEC& debug1, VEC& debug2, VEC& debug3,
     Dune::ParameterTree& ptree)
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
    //     a) solving the local geneo eigenproblem and picking all eigenvectors with eigenvalue smaller than given threshold
    //     b) multiplying these eigenvectors with the partition of unity
    //     c) nulling out dofs on the global dirchlet boundary
    // this way, a prolongation to the global Dirichlet boundary will produce a zero contribution
    //
    // this is the only step that was replaced to turn TwoLevelSchwarz into TwoLevelGenEO!
    //

    using Dune::PDELab::Backend::native;

    //auto n_coarse_basis_vectors_per_subdomain = input_coarse_vectors.size(); // here the number of coarse vectors in each subdomain is the same
    for (unsigned p=0; p<subdomains; p++)
    {
       // set up matrices

      // A is the neumann matrix. It will be the left-hand side of our eigenproblem. (This is redundant but writing 'A' is nicer.)
      auto A = (*matrices)[p]; 

      // B will be the right-hand side of our eigenproblem. It starts as a physical copy of the neumann matrix.
      auto B = *(*matrices)[p];

      // Now multiply B from left and right with the partition of unity.
      for (auto row_iter = B.begin(); row_iter != B.end(); row_iter++) 
      {
        for (auto col_iter = row_iter->begin(); col_iter != row_iter->end(); col_iter++) 
        {
          *col_iter *= partition_of_unity[p][row_iter.index()] * partition_of_unity[p][col_iter.index()];
        }
      }

      // // write matrices to file for debugging
      // if (p == 7)
      // {
      //   Dune::storeMatrixMarket(*A,"process_" + std::to_string(p) + "_A.mtx");
      //   Dune::storeMatrixMarket(B,"process_" + std::to_string(p) + "_B.mtx");
      // }

      // Setup ArpackGenEO wrapper, solving a generalized eigenproblem with lhs A.
      ArpackGeneo::ArPackPlusPlus_Algorithms<Dune::PDELab::Backend::Native<MAT>, VEC> arpack(*A);

      // set some parameters. Some or all of these should be input parameters.
      double tolerance = ptree.get<double>("geneo.arpack_tolerance",1e-6); // tolerance for arpack algorithm
      int number_of_eigenvalues = ptree.get<int>("geneo.number_of_eigenvectors_arpack",10); // how many eigenvalues does arpack compute
      double threshold = ptree.get<double>("geneo.threshold",-1.0); // eigenvalues below this are considered interesting. Their eigenvectors are added to coarse_basis_vectors. If <0, then use fixed number given below.
      int number_of_basisvectors = ptree.get<int>("geneo.number_of_basisvectors",3); // number of eigenvectors added to coarse_basis_vectors if threshold < 0
      std::cout << p << ": Number of arpack eigenpairs = " << number_of_eigenvalues <<" Number of basisvectors used = " << number_of_basisvectors << std::endl;
      double shift = 0.001; // shift used in arpack algorihm
      int verbose = 1;

      // make sure to allocate the right amount of space to our coarse_basis_vectors on subdomain [p]
      auto n = (*local_to_global)[p].size(); // size of the coarse vectors in subdomain p (in blocks)
      VEC vec(n); // make a template
      coarse_basis_vectors[p].resize(number_of_eigenvalues,vec); // create the subdomain istl vectors; now coarse_basis_vectors[p][k] is an ISTL vector on subdomain p        

      // make vectors to store eigenvalues. Eigenvectors will be written directly into coarse_basis_vectors.
      std::vector<double> eigenvalues(number_of_eigenvalues,0.0);

      // solve GenEO eigenproblem
      arpack.computeGenNonSymMinMagnitude(B, tolerance, coarse_basis_vectors[p], eigenvalues, shift); 


      // // Count eigenvectors below threshold
      // if (threshold >= 0) {
      //   for (int i = 0; i < number_of_eigenvalues; i++) {
      //     if (eigenvalues[i] > threshold) {
      //       number_of_basisvectors = i;
      //       break;
      //     }
      //   }
      //   if (verbose > 0)
      //     std::cout << "Process " << p << " picked " << number_of_basisvectors << " eigenvectors" << std::endl;
      //   if (number_of_basisvectors == -1)
      //     DUNE_THROW(Dune::Exception,"No eigenvalue above threshold - not enough eigenvalues computed!");
      // } else {
      //   number_of_basisvectors = number_of_eigenvalues;
      // }

      // resize coarse_basis_vectors to accomodate only number_of_basisvectors basis vectors, deleting the rest
      coarse_basis_vectors[p].resize(number_of_basisvectors,vec); // create the subdomain istl vectors; now coarse_basis_vectors[p][k] is an ISTL vector on subdomain p        

    	// modify coarse_basis_vectors
    	for (size_t k=0; k<number_of_basisvectors; k++)
      {
        field_type norm = 0.0;
        for (size_t i=0; i<n; i++)
        {
          norm += coarse_basis_vectors[p][k][i].two_norm2();
        }
        norm = sqrt(norm);
    	  for (size_t i=0; i<n; i++)
  	    {
          coarse_basis_vectors[p][k][i] *= 1.0/norm;
  	      auto iglobal = (*local_to_global)[p][i];  
  	      coarse_basis_vectors[p][k][i] *= partition_of_unity[p][i]; // multiply with partition of unity, puts zero on subdomain boundary
  	      for (size_t j=0; j<blocksize; j++)
  		      coarse_basis_vectors[p][k][i][j] *= (*global_dirichlet_dofs)[iglobal][j]; // zero out dof on global Dirichlet boundary	      
  	    }
      }


      // if (p ==7) // write some output, but only for one subdomain
      // {
      //   // write eigenvalues into output for debugging
      //   for (size_t i = 0; i < eigenvalues.size(); i++)
      //     {
      //       std::cout << "Eigenvalue number " << i << " = " << eigenvalues[i] << std::endl;
      //     }
      
      //   // As a test for the neumann matrix A, set up a constant function (an istl vector containing only ones). 
      //   // This should be an eigenfunction of $A$ with eigenvalue 0.

      //   // vector representing constant function
      //   VEC constant_vector(n);
      //   for (size_t i=0; i<n; i++)
      //   {
      //     auto iglobal = (*local_to_global)[p][i];
      //     constant_vector[i] = (*input_coarse_vectors[0])[iglobal]; // write 0-th coarse input vector into 'constant_vector'
      //   }

      //   // set up a vector to store A*constant_vector, where constant one is a vector containing ones.
      //   VEC A_times_constant(n);
      //   A_times_constant = 0; //set to 0

      //   // calculate A_times_constant = A * constant_vector 
      //   auto block = constant_vector[0]; // temporary storage for a block, used in the loop below
      //   for (auto row_iter = A->begin(); row_iter != A->end(); row_iter++) 
      //   {
      //     for (auto col_iter = row_iter->begin(); col_iter != row_iter->end(); col_iter++) {
      //       col_iter->mv(constant_vector[col_iter.index()], block); // multiply matrix block with vector block
      //       A_times_constant[row_iter.index()] += block;
      //     }
      //   }

      //   // // use constant vector times oartition of unity as first coarse basis vector (to have it as output.)
      //   // coarse_basis_vectors[p][0] = A_times_constant;
      //   // coarse_basis_vectors[p][1] = constant_vector;

      //   // for (size_t i=0; i<n; i++)
      //   // {
      //   //   auto iglobal = (*local_to_global)[p][i];  
      //   //   coarse_basis_vectors[p][0][i] *= partition_of_unity[p][i]; // multiply with partition of unity, puts zero on subdomain boundary
      //   //   for (size_t j=0; j<blocksize; j++)
      //   //     coarse_basis_vectors[p][0][i][j] *= (*global_dirichlet_dofs)[iglobal][j]; // zero out dof on global Dirichlet boundary        
      //   // }

      //   // calculate norm of A_times_constant
      //   double norm = 0.0;
      //   for (size_t i = 0; i < n; i++)
      //   {
      //     norm += A_times_constant[i]*A_times_constant[i];
      //   }

      //   std::cout << "||A * constant||^2 = " << norm << std::endl;
      // } // end if

      // overwrite them again
      // for (size_t k=3; k<number_of_basisvectors; k++)
      //   for (size_t i=0; i<n; i++)
      //   {
      //     auto iglobal = (*local_to_global)[p][i];
      //     coarse_basis_vectors[p][k][i] = (*input_coarse_vectors[k])[iglobal]; // pick out
      //     coarse_basis_vectors[p][k][i] *= partition_of_unity[p][i]; // multiply with partition of unity, puts zero on subdomain boundary
      //     for (size_t j=0; j<blocksize; j++)
      //       coarse_basis_vectors[p][k][i][j] *= (*global_dirichlet_dofs)[iglobal][j]; // zero out dof on global Dirichlet boundary        
      //   }
      //coarse_basis_vectors[p].resize(input_coarse_vectors.size());

    } // end for loop over subdomains


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
    int show1=0, show2=0;
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
	
	// now do the triple matrix products
	VEC v1(Asubdomain.N()); // two vectors on the subdomain to compute scalar products
	VEC v2(Asubdomain.N());
	for (size_t krow=0; krow<coarse_basis_vectors[prow].size(); krow++) // loop over all coarse grid functions in this subdomain
	  {
	    // row number in coarse system
	    auto icoarse = coarse_vector_local_to_global[prow][krow];
	    
	    // matrix vector product
	    Asubdomain.mv(coarse_basis_vectors[prow][krow],v1);

	    // debugging vectors
	    if (prow==show1 && krow==2)
	      {
		debug1 = 0.0;
		debug2 = 0.0;
		for (size_t i=0; i<Asubdomain.N(); i++)
		  {
		    debug1[(*local_to_global)[prow][i]] = coarse_basis_vectors[prow][krow][i];
		    debug2[(*local_to_global)[prow][i]] = v1[i];
		  }
	      }
	      

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

		  // debug output
		  if (prow==show1 && pcol==show2 && krow==2 && kcol==0)
		    {
		      debug3 = 0.0;
		      for (size_t i=0; i<Asubdomain.N(); i++)
			debug3[(*local_to_global)[prow][i]] = v2[i];
		    }

		  // now we have a matrix entry
		  (*pcmat)[icoarse][jcoarse] = v2*v1;
		}
	  }
      }
    std::cout << "computed entries of coarse system" << std::endl;
    // Dune::printmatrix(std::cout, *pcmat, "coarse level matrix", "");

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

  // get a prolongated coarse basis vector (to write to vtk for debugging)
  // p = subdomain number
  // k = number of basisvector from this subdomain
  VEC prolongated_basis_vector(int p, int k)
  {
    auto n_local = (*local_to_global)[p].size(); // size of local vector
    auto n_global = (*global_to_local).size(); // size of global vector

    VEC basis_vector(n_global); //initialize basis_vector with global size

    // fill basis_vector with zeros
    for (size_t i = 0; i<basis_vector.size(); i++)
    {
      for (size_t j = 0; j<basis_vector[i].size(); j++)
      {
        basis_vector[i][j] = 0.0;
      }
    }

    // copy coarse_basis_vector[p][k] into basis_vector at the appropriate indices
    for (size_t i = 0; i < n_local; i++)
    {
      basis_vector[(*local_to_global)[p][i]] = coarse_basis_vectors[p][k][i];
    }

    return basis_vector;
  }

};

#endif // HAVE_SUITESPARSE_UMFPACK || DOXYGEN

#endif // Udune_ftworth_HH
