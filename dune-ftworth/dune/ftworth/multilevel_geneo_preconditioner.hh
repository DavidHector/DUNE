
#ifndef Udune_ftworth_multilevel_geneo_preconditioner_HH
#define Udune_ftworth_multilevel_geneo_preconditioner_HH

#include "partitioner.hh"
#include "subdomainutilities.hh"
#include "schwarz_preconditioners.hh"
#include "arpack_geneo_wrapper.hh"

#define USE_CHOLMOD 1

#if (HAVE_SUITESPARSE_UMFPACK && HAVE_EIGEN) || DOXYGEN

#include <dune/istl/umfpack.hh>
#include <dune/istl/cholmod.hh>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

template<typename field_type>
std::vector<std::vector<field_type>>
get_partition_of_unity_standard (const std::vector<std::vector<size_t>>& local_to_global,
                                 const std::vector<std::map<size_t,size_t>>& global_to_local,
                                 const std::vector<std::vector<bool>>& boundary_dofs)
{
  // number of subdomains we have
  auto subdomains = local_to_global.size();

  // make a global vector holding the number of subdomains where dof is NOT on the boundary
  std::vector<int> k(global_to_local.size(),0);
  for (unsigned p=0; p<subdomains; p++)
    for (unsigned i=0; i<local_to_global[p].size(); i++)
      if (!boundary_dofs[p][i])
        k[local_to_global[p][i]] += 1;

  // now compute partition of unity
  std::vector<std::vector<field_type>> pu(subdomains);
  for (unsigned p=0; p<subdomains; p++)
    {
      pu[p].resize(local_to_global[p].size());
      for (unsigned i=0; i<local_to_global[p].size(); i++)
        pu[p][i] = (boundary_dofs[p][i]) ? 0.0 : (((field_type)1.0)/k[local_to_global[p][i]]);
    }

  // return
  return pu;
}

template<typename MAT>
std::vector<std::vector<typename MAT::field_type>>
get_partition_of_unity_distance (const std::vector<std::shared_ptr<MAT>>& matrices,
                                 const std::vector<std::vector<size_t>>& local_to_global,
                                 const std::vector<std::map<size_t,size_t>>& global_to_local,
                                 const std::vector<std::vector<bool>>& boundary_dofs)
{
  // types
  using field_type = typename MAT::field_type;

  // number of subdomains we have
  auto subdomains = local_to_global.size();

  // compute local distance from boundary in each subdomain
  const int maxdistance = 1<<24;
  std::vector<std::vector<int>>localdistance(subdomains);
  for (unsigned p=0; p<subdomains; p++)
    {
      // dof blocks in subdomain p
      auto n = local_to_global[p].size();

      // allocate subdomain vector and initialize with a distance never reached
      localdistance[p].resize(n,maxdistance);

      // initialize boundary degrees of freedom and interior degrees of freedom
      for (unsigned i=0; i<n; i++)
        {
          if (boundary_dofs[p][i]) localdistance[p][i] = 0; // boundary case
          if ( global_to_local[local_to_global[p][i]].size()==1 ) localdistance[p][i] = -1; // interior dof case
        }
      // non boundary and non interior dofs are now at maxdistance

      // get the subdomain matrix
      MAT& A = *matrices[p];

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
  std::vector<int> globaldistance(global_to_local.size(),0);
  for (unsigned p=0; p<subdomains; p++)
    for (unsigned i=0; i<local_to_global[p].size(); i++)
      if (localdistance[p][i]>=0)
        globaldistance[local_to_global[p][i]] += localdistance[p][i];

  // now compute partition of unity
  std::vector<std::vector<field_type>> pu(subdomains);
  for (unsigned p=0; p<subdomains; p++)
    {
      pu[p].resize(local_to_global[p].size());
      for (unsigned i=0; i<local_to_global[p].size(); i++)
        {
          if (localdistance[p][i]>=0 && globaldistance[local_to_global[p][i]]==0)
            std::cout << "compute PU zero global distance! p=" << p << " i=" << i << " g=" << local_to_global[p][i] << " local distance=" << localdistance[p][i] << " global distance=" << globaldistance[local_to_global[p][i]] << std::endl;
          if (localdistance[p][i]<0)
            pu[p][i] = 1.0;
          else
            pu[p][i] = ((field_type)localdistance[p][i])/globaldistance[local_to_global[p][i]];
        }
    }

  // return
  return pu;
}


template<typename field_type>
std::vector<std::vector<field_type>>
get_subdomain_boundary_mask (const std::vector<std::vector<bool>>& boundary_dofs)
{
  auto subdomains = boundary_dofs.size();
  std::vector<std::vector<field_type>> boundary_mask(subdomains);
  for (size_t p=0; p<subdomains; p++)
    {
      auto n = boundary_dofs[p].size(); // dof blocks in subdomain p
      boundary_mask[p].resize(n);
      for (size_t i=0; i<n; i++)
        boundary_mask[p][i] = (boundary_dofs[p][i]) ? 0.0 : 1.0;
    }
  return boundary_mask;
}


/** \brief Multi-level GenEO preconditioner
 *
 * \author Peter Bastian
 *
 * \tparam MAT a "Dune::BCRSMatrix<Dune::FieldMatrix>" type
 * \tparam VEC a "Dune::BlockVector<Dune::FieldVector>" type
 *
 */
template<typename MAT, typename VEC>
class MultiLevelGenEO : public Dune::Preconditioner<VEC,VEC>
{
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

private:
  // types use din this class
  using VectorBlockType = typename VEC::block_type;
  using CBLOCKMAT = Dune::FieldMatrix<typename VEC::field_type,1,1>;
  using CBLOCKVEC = Dune::FieldVector<typename VEC::field_type,1>;
  using DenseMat = Dune::DynamicMatrix<field_type>; // type for coarse snippet matrices
  using CMAT = Dune::BCRSMatrix<CBLOCKMAT>; // type for coarse level matrix
  using CVEC = Dune::BlockVector<CBLOCKVEC>; // type for coarse level vector
  using CHOLMODSOLVER = Dune::Cholmod<VEC>;
  using CCHOLMODSOLVER = Dune::Cholmod<CVEC>;
  using UMFPACKSOLVER = Dune::UMFPack<MAT>;
  using CUMFPACKSOLVER = Dune::UMFPack<CMAT>;

  enum { blocksize = VectorBlockType::dimension };

  // the local data of this class

  // input to constructor
  MPI_Comm comm;
  std::shared_ptr<MAT> pAglobal;
  std::shared_ptr<VEC> pglobal_dirichlet_mask;
  std::shared_ptr<std::vector<std::shared_ptr<MAT>>> pvolume_snippet_matrices;
  std::shared_ptr<std::vector<std::shared_ptr<MAT>>> pskeleton_snippet_matrices;
  std::shared_ptr<DomainDecompositionDOFMapper> pdddm;
  std::vector<size_t> nsubdomains;

  // local data
  size_t Number_Of_Levels; // how many levels do we have. Note 0 is the *finest* level and Number_Of_Levels-1 is the *coarsest* level!
  std::vector<std::shared_ptr<DomainDecomposition>> ddhierarchy; // hierarchy of domain decomposition information
  std::vector<std::shared_ptr<DomainDecompositionDOFMapper>> dddmhierarchy; // hierarchy of dof information on the levels

  // system information
  std::vector<std::vector<std::vector<bool>>> subdomain_boundary_dofs; // subdomain_boundary_dofs[l][p][i] is true if dof i in subdomain p on level l is on the boundary
  std::vector<std::vector<std::vector<field_type>>> subdomain_boundary_mask; // subdomain_boundary_dofs converted to a mask
  std::vector<std::vector<std::vector<field_type>>> partition_of_unity; // partition_of_unity[l][p][i] is in [0,1] for dof i in subdomain p on level l
  std::vector<std::vector<VEC>> fine_subdomain_basis_vectors; //  basis vectors computed in each subdomain on the finest level
  std::vector<std::vector<std::vector<CVEC>>> coarse_subdomain_basis_vectors; //  basis vectors computed in each subdomain on the coarser level coarse_subdomain_basis_vectors[l][p][k]
  std::vector<std::vector<std::vector<CVEC>>> coarse_kerAB_vectors; //  basis vectors computed in each subdomain on the coarser level coarse_subdomain_basis_vectors[l][p][k]
  std::vector<std::shared_ptr<CMAT>> global_coarse_matrix; // the global systems on the coarse levels
  std::vector<std::vector<std::shared_ptr<DenseMat>>> projected_volume_snippet_matrices; // volume snippet matrices projected to subdomain basis vectors
  std::vector<std::vector<std::shared_ptr<DenseMat>>> projected_skeleton_snippet_matrices; // skeleton snippet matrices projected to subdomain basis vectors
  std::vector<std::vector<std::shared_ptr<CMAT>>> volume_snippet_matrices; // volume snippet matrices on coarser levels; from these the neumann matrices are assembled
  std::vector<std::vector<std::shared_ptr<CMAT>>> skeleton_snippet_matrices; // skeleton snippet matrices on coarser levels; from these the neumann matrices are assembled

  // subdomain level solves
#ifdef USE_CHOLMOD
  std::vector<std::shared_ptr<CHOLMODSOLVER>> finesubdomainsolvers; // the LU decomposition of the subdomain problems on the finest level
  std::vector<std::shared_ptr<CCHOLMODSOLVER>> globalcoarsesolvers; // the LU decomposition of the coarse level problems
  std::vector<std::vector<std::shared_ptr<CCHOLMODSOLVER>>> coarsesubdomainsolvers; // the subdomain solvers on coarse levels
#else
  std::vector<std::shared_ptr<UMFPACKSOLVER>> finesubdomainsolvers; // the LU decomposition of the subdomain problems on the finest level
  std::vector<std::shared_ptr<CUMFPACKSOLVER>> globalcoarsesolvers; // the LU decomposition of the coarse level problems
  std::vector<std::vector<std::shared_ptr<CUMFPACKSOLVER>>> coarsesubdomainsolvers; // the subdomain solvers on coarse levels
#endif
  std::vector<std::shared_ptr<CVEC>> dcoarse; //  restricted defect
  std::vector<std::shared_ptr<CVEC>> vcoarse; //  prolongated corrections

  // timing
  std::vector<std::vector<double>> setup_time_procs; // time spent on level l in proc p
  std::vector<std::vector<double>> setup_time_setup_gevp; // time spent on level l in proc p
  std::vector<std::vector<double>> setup_time_QR; // time spent on level l in proc p
  std::vector<std::vector<double>> setup_time_projection; // time spent on level l in proc p
  std::vector<std::vector<double>> setup_time_gevp; // time spent on level l in proc p
  std::vector<std::vector<double>> setup_time_umfpack; // time spent on level l in proc p
  std::vector<double> setup_time_parallelizable; // times now sequential but probably parallelizable on each level
  double setup_time_sequential; // really sequential
  double setup_time_total;

  // cycle params
  std::string cycle; // controls cycle form: "additive", "multiplicativeras"

public:
  //! \brief Constructor. It does all the work of setting up the components
  /**
   * \param comm_                        an MPI communicator (only needed for parmetis)
   * \param pAglobal_                    shared_ptr to the global system matrix
   * \param pglobal_dirichlet_mask_      global dof vector with 0 in dirichlet dofs and 1 else
   * \param pvolume_snippet_matrices_    assembled volume snippet matrices on finest level
   * \param pskeleton_snippet_matrices_  assembled skeleton snippet matrices on finest level
   * \param pdddm_                       dof mapper and domain decomposition information on finest level
   * \param nsubdomains                  number of subdomains on the different levels
   * \param coarse_overlap_              overlap to be used when building coarse level domain decomposition
   * \param pum                          partition of unity method to be used on finest level
   */
  MultiLevelGenEO (MPI_Comm comm_, // an MPI communicator (only needed for parmetis)
                   std::shared_ptr<MAT> pAglobal_, // the global stiffness matrix
                   std::shared_ptr<VEC> pglobal_dirichlet_mask_, // global dof vector with 0 in dirichlet dofs and 1 else
                   std::shared_ptr<std::vector<std::shared_ptr<MAT>>> pvolume_snippet_matrices_, // assembled volume snippet matrices on finest level
                   std::shared_ptr<std::vector<std::shared_ptr<MAT>>> pskeleton_snippet_matrices_, // assembled skeleton snippet matrices on finest level
                   std::shared_ptr<DomainDecompositionDOFMapper> pdddm_, // dof mapper and domain decomposition information on finest level
                   const std::vector<size_t>& nsubdomains_, // desired subdomain numbers on all levels
                   std::string pum_, // pu to be used
                   std::string fineGEVPrhs_, // variant of B: use "pu", "1-pu" with "pu" being the default
                   std::string coarseGEVPrhs_, // variant of B: use "pu", "1-pu" with "pu" being the default
                   size_t coarse_overlap_, // overlap to be used when building coarse level domain decomposition
                   std::string coarseeigensolver_, // can take either "eigen" or "arpack" with eigen as the default
                   field_type arpack_tolerance_, // threshold for iterative eigensolver ARPACK, use 0.0 for internal default
                   size_t n_eigenvectors_fine_computed_, size_t n_eigenvectors_fine_used_,
                   size_t n_eigenvectors_coarse_computed_, size_t n_eigenvectors_coarse_used_,
                   field_type eigenvalue_fine_threshold_, field_type eigenvalue_coarse_threshold_,// <=0 means fixed size, >0 means take only evs with eiganvale<=threshold
                   field_type abs_zero_ker_,
                   field_type regularization_ker_,
                   bool merge_disconnected_,
                   std::string cycle_,
                   size_t verbose_ // 0=no output during creation, 1=minimal information
                   )
    : comm(comm_),
      pAglobal(pAglobal_),
      pglobal_dirichlet_mask(pglobal_dirichlet_mask_),
      pvolume_snippet_matrices(pvolume_snippet_matrices_),
      pskeleton_snippet_matrices(pskeleton_snippet_matrices_),
      pdddm(pdddm_),
      nsubdomains(nsubdomains_),
      cycle(cycle_)
  {
    Dune::Timer timer_total, timer, timer2; // two timers, one for total time, one for general use
    timer_total.reset();

    // some analysis of the input
    if (verbose_>0) std::cout << "starting level " << 0 << std::endl;
    size_t totaldofs=blocksize*pglobal_dirichlet_mask->N(), constraineddofs=0;
    for (size_t i=0; i<pglobal_dirichlet_mask->N(); i++)
      for (size_t j=0; j<blocksize; j++)
        if ( (*pglobal_dirichlet_mask)[i][j] == 0.0 )
          constraineddofs++;
    std::cout << "ML Geneo blocksize=" << blocksize << " constraineddofs=" << constraineddofs << " out of " << totaldofs << std::endl;

    // how many levels do we have. Note 0 is the *finest* level and Number_Of_Levels-1 is the *coarsest* level!
    Number_Of_Levels = nsubdomains.size();

    // check input
    if (Number_Of_Levels<2)
      DUNE_THROW(Dune::Exception, "number of levels must be at least 2.");
    if (nsubdomains[Number_Of_Levels-1]!=1)
      DUNE_THROW(Dune::Exception, "number of subdomains on coarsest level must be 1.");

    //====================
    // We start on the finest level
    //====================

    // allocate time slots
    setup_time_procs.resize(Number_Of_Levels);
    setup_time_setup_gevp.resize(Number_Of_Levels);
    setup_time_QR.resize(Number_Of_Levels);
    setup_time_projection.resize(Number_Of_Levels);
    setup_time_gevp.resize(Number_Of_Levels);
    setup_time_umfpack.resize(Number_Of_Levels);
    setup_time_procs[0].resize(nsubdomains[0],0.0);
    setup_time_umfpack[0].resize(nsubdomains[0],0.0);
    setup_time_parallelizable.resize(Number_Of_Levels,0.0);
    setup_time_sequential = 0.0;
    setup_time_total = 0.0;

    // initialize multilevel domain decomposition and associated dof information
    ddhierarchy.resize(Number_Of_Levels);
    dddmhierarchy.resize(Number_Of_Levels);
    ddhierarchy[0] = pdddm->pdd;
    ddhierarchy[0]->print_info(verbose_);
    dddmhierarchy[0] = pdddm;
    dddmhierarchy[0]->print_info(verbose_);

    // for the solver application and the partition of unity we need the subdomain boundary dof information
    subdomain_boundary_dofs.resize(Number_Of_Levels);
    subdomain_boundary_dofs[0] = get_subdomain_boundary_dofs(*pAglobal,pdddm->subdomainnumber_local_to_global,pdddm->subdomainnumber_global_to_local);
    for (size_t p=0; p<nsubdomains[0]; ++p)
      {
        size_t sdboundarydofs=0;
        for (size_t i=0; i<subdomain_boundary_dofs[0][p].size(); i++)
          if (subdomain_boundary_dofs[0][p][i])
            sdboundarydofs += blocksize;
        std::cout << "ML Geneo p=" << p << " boundarydofs=" << sdboundarydofs << std::endl;
      }

    // convert boundary dof information to a mask
    subdomain_boundary_mask.resize(Number_Of_Levels);
    subdomain_boundary_mask[0] = get_subdomain_boundary_mask<field_type>(subdomain_boundary_dofs[0]);

    // build partition of unity on finest level
    // we can handle only the standard one here; otherwise we would have to assemble all subdomain matrices at once
    // this fits also better to the coarser levels which can only used the standard pu
    timer.reset();
    partition_of_unity.resize(Number_Of_Levels);
    if (pum_=="distance")
      {
        // we need all matrices, this is costly!
        auto matrices = set_up_local_matrices(*pAglobal,pdddm->subdomainnumber_local_to_global,pdddm->subdomainnumber_global_to_local);
        partition_of_unity[0] = get_partition_of_unity_distance(matrices,pdddm->subdomainnumber_local_to_global,pdddm->subdomainnumber_global_to_local,subdomain_boundary_dofs[0]);
      }
    else
      partition_of_unity[0] = get_partition_of_unity_standard<field_type>(pdddm->subdomainnumber_local_to_global,pdddm->subdomainnumber_global_to_local,subdomain_boundary_dofs[0]);
    setup_time_parallelizable[0] += timer.elapsed();

    // now lets work on the finest level: assemble subdomain matrices, etc
    finesubdomainsolvers.resize(nsubdomains[0]); // we will store the UMFPack solver objects here
    fine_subdomain_basis_vectors.resize(nsubdomains[0]);

    // it suffices to assemble subdomain problems one by one!
    std::vector<bool> disconnected_subdomains(nsubdomains[0],false); // true when subdomain is not connected
    for (size_t p=0; p<nsubdomains[0]; ++p)
      {
        timer.reset();

        // set up sparse subdomain matrix
        std::shared_ptr<MAT> pA = set_up_local_matrix(*pAglobal,p,pdddm->subdomainnumber_local_to_global,pdddm->subdomainnumber_global_to_local);
        *pA = 0.0; // clear subdomain matrix

        // assemble snippets to subdomain matrix
        // This results in Neumann subdomain matrix without any boundary conditions!
        assemble_snippets_to_subdomain_matrix(p,
                                              pdddm->pdd->subdomain_to_volumesnippetnumber,
                                              pdddm->volumesnippetnumber_local_to_global,
                                              pdddm->subdomainnumber_global_to_local,
                                              *pvolume_snippet_matrices,
                                              pA);
        if (pskeleton_snippet_matrices!=nullptr)
          assemble_snippets_to_subdomain_matrix(p,
                                                pdddm->pdd->subdomain_to_interior_skeletonsnippetnumber,
                                                pdddm->skeletonsnippetnumber_local_to_global,
                                                pdddm->subdomainnumber_global_to_local,
                                                *pskeleton_snippet_matrices,
                                                pA);


        // assemble global dirichlet constraints in subdomain matrix
        bool floating = assemble_dirichlet_constraints_in_subdomain(pdddm->subdomainnumber_local_to_global[p],*pglobal_dirichlet_mask,pA);

        // now at this point we have the NEUMANN subdomain matrix (with global Dirichlet conditions enforced in non-floating subdomains)
        // now we may set up and solve the eigenvalue problems
        MAT& A = *pA; // introduce reference for ease of writing
        make_symmetric(A); // out of paranoia

        // for (auto row_iter = pA->begin(); row_iter != pA->end(); ++row_iter)
        //   for (size_t i=0; i<blocksize; i++)
        //     {
        //       auto double maxvalue = 0.0;
        //       for (auto col_iter = row_iter->begin(); col_iter != row_iter->end(); ++col_iter)
        //      for (size_t j=0; j<blocksize; j++)
        //        maxvalue = std::max(maxvalue,std::abs((*col_iter)[i][j]));

        //       for (auto col_iter = row_iter->begin(); col_iter != row_iter->end(); ++col_iter)
        //      for (size_t j=0; j<blocksize; j++)
        //        {
        //          auto rowglobal = pdddm->subdomainnumber_local_to_global[p][row_iter.index()];
        //          auto colglobal = pdddm->subdomainnumber_local_to_global[p][col_iter.index()];
        //          auto aijlocal = (*col_iter)[i][j];
        //          auto aijglobal = (*pAglobal)[rowglobal][colglobal][i][j];
        //          if ( (*pglobal_dirichlet_mask)[rowglobal][i]!=0.0 && (*pglobal_dirichlet_mask)[colglobal][j]!=0.0
        //               && !subdomain_boundary_dofs[0][p][row_iter.index()] && !subdomain_boundary_dofs[0][p][col_iter.index()] )
        //            {
        //              if (maxvalue<1e-13 && std::abs(aijlocal-aijglobal)>1e-13)
        //                std::cout << "p=" << p << " rowlocal=" << row_iter.index() << " collocal=" << col_iter.index()
        //                          << " rowglobal=" << rowglobal << " colglobal=" << colglobal << " i=" << i << " j=" << j
        //                          << " aijl=" << aijlocal << " aijg=" << aijglobal << " absdiff=" << std::abs(aijlocal-aijglobal) << std::endl;
        //              if (maxvalue>=1e-13 && std::abs(aijlocal-aijglobal)/maxvalue>1e-13)
        //                std::cout << "p=" << p << " rowlocal=" << row_iter.index() << " collocal=" << col_iter.index()
        //                          << " rowglobal=" << rowglobal << " colglobal=" << colglobal << " i=" << i << " j=" << j
        //                          << " aijl=" << aijlocal << " aijg=" << aijglobal << " reldiff=" << std::abs(aijlocal-aijglobal)/maxvalue << std::endl;
        //            }
        //        }
        //     }

        // B will be the right-hand side of our eigenproblem. It starts as a physical copy of the neumann matrix.
        MAT B(A);

        // Now multiply B from left and right with the partition of unity.
        std::string basename("ftworth");
        if (fineGEVPrhs_=="geneo")
          {
            B = 0.0; // clear B
            assemble_snippets_to_subdomain_matrix_geneo_volume(p,
                                                               pdddm->pdd->subdomain_to_volumesnippetnumber,
                                                               pdddm->pdd->number_to_volumesnippet,
                                                               pdddm->volumesnippetnumber_local_to_global,
                                                               pdddm->subdomainnumber_global_to_local,
                                                               *pvolume_snippet_matrices,
                                                               B); // assemble without interior
            if (pskeleton_snippet_matrices!=nullptr)
              assemble_snippets_to_subdomain_matrix_geneo_skeleton(p,
                                                                   pdddm->pdd->subdomain_to_interior_skeletonsnippetnumber,
                                                                   pdddm->pdd->number_to_volumesnippet,
                                                                   pdddm->pdd->number_to_skeletonsnippet,
                                                                   pdddm->skeletonsnippetnumber_local_to_global,
                                                                   pdddm->subdomainnumber_global_to_local,
                                                                   *pskeleton_snippet_matrices,
                                                                   B);
            assemble_dirichlet_constraints_in_subdomain(pdddm->subdomainnumber_local_to_global[p],*pglobal_dirichlet_mask,B);

            // std::ostringstream sBovlp;
            // sBovlp << basename << "_Bovlp_" << 1 << "_" << p << ".txt";
            // std::cout << "writing matrix file " << sBovlp.str() << std::endl;
            // writeMatrixToMatlab(B,sBovlp.str(),16);

            for (auto row_iter = B.begin(); row_iter != B.end(); ++row_iter)
              for (auto col_iter = row_iter->begin(); col_iter != row_iter->end(); ++col_iter)
                *col_iter *= partition_of_unity[0][p][row_iter.index()] * partition_of_unity[0][p][col_iter.index()];
          }
        else if (fineGEVPrhs_=="1-pu")
          {
            for (auto row_iter = B.begin(); row_iter != B.end(); ++row_iter)
              for (auto col_iter = row_iter->begin(); col_iter != row_iter->end(); ++col_iter)
                *col_iter *= (1.0-partition_of_unity[0][p][row_iter.index()]) * (1.0-partition_of_unity[0][p][col_iter.index()]);
          }
        else
          {
            for (auto row_iter = B.begin(); row_iter != B.end(); ++row_iter)
              for (auto col_iter = row_iter->begin(); col_iter != row_iter->end(); ++col_iter)
                *col_iter *= partition_of_unity[0][p][row_iter.index()] * partition_of_unity[0][p][col_iter.index()];
          }

        // writing matrices
        // std::ostringstream sA;
        // sA << basename << "_A_" << 1 << "_" << p << ".txt";
        // std::cout << "writing matrix file " << sA.str() << std::endl;
        // writeMatrixToMatlab(A,sA.str(),16);
        // std::ostringstream sB;
        // sB << basename << "_B_" << 1 << "_" << p << ".txt";
        // std::cout << "writing matrix file " << sB.str() << std::endl;
        // writeMatrixToMatlab(B,sB.str(),16);
        // MAT D(A);
        // D = 0.0;
        // for (auto row_iter = D.begin(); row_iter != D.end(); ++row_iter)
        //   for (size_t i=0; i<blocksize; i++)
        //     D[row_iter.index()][row_iter.index()][i][i] = partition_of_unity[0][p][row_iter.index()];
        // std::ostringstream sD;
        // sD << basename << "_D_" << 1 << "_" << p << ".txt";
        // std::cout << "writing matrix file " << sD.str() << std::endl;
        // writeMatrixToMatlab(D,sD.str(),16);

        // Setup ArpackGenEO wrapper, solving a generalized eigenproblem with lhs A.
        ArpackMLGeneo::ArPackPlusPlus_Algorithms<MAT,VEC> arpack(A);

        // set some parameters. Some or all of these should be input parameters.
        double tolerance = arpack_tolerance_; // tolerance for arpack algorithm
        double shift = 0.001; // shift used in arpack algorihm

        // make seperate vectors for arpack
        auto n = pdddm->subdomainnumber_local_to_global[p].size(); // size of vectors in subdomain p
        VEC vec(n); // make a template
        vec = 0.0;
        std::vector<VEC> eigenvectors(n_eigenvectors_fine_computed_,vec); // create vector of ISTL vectors for Arpack to compute

        // make vectors to store eigenvalues. Eigenvectors will be written directly into coarse_basis_vectors.
        std::vector<double> eigenvalues(n_eigenvectors_fine_computed_,0.0);

        // solve GenEO eigenproblems
        //if (verbose_>0)
        std::cout << "subdomain " << p << ": solve eigenvalue problem, floating=" << floating << std::endl;

        //arpack.computeStdNonSymMinMagnitude(B,tolerance,eigenvectors,eigenvalues,shift);
        //arpack.computeGenSymShiftInvertMinMagnitude(B,tolerance,eigenvectors,eigenvalues,shift);
        if (eigenvalue_fine_threshold_<=0)
          arpack.computeGenSymShiftInvertMinMagnitude(B,tolerance,eigenvectors,eigenvalues,shift);
        else
          arpack.computeGenSymShiftInvertMinMagnitudeAdaptive(B,tolerance,eigenvectors,eigenvalues,shift,eigenvalue_fine_threshold_,n_eigenvectors_fine_used_);
        //if (verbose_>-1)
          for (size_t i=0; i<eigenvalues.size(); ++i)
            std::cout << "  lambda_" << i<< " = " << eigenvalues[i] << std::endl;

        // analyse eigenvalues
        if (eigenvalues[0]<-1e-4)
          {
            std::cout << "negative eigenvalue(s) found in subdomain " << p << " on level " << 0 << std::endl;
            for (size_t i=0; i<eigenvalues.size(); ++i)
              std::cout << "  lambda_" << i<< " = " << eigenvalues[i] << std::endl;
            DUNE_THROW(Dune::Exception, "at least one negative eigenvalue detected");
          }
        size_t number_of_zero_eigenvalues = 0;
        // for (size_t i=0; i<eigenvalues.size(); ++i)
        //   if (std::abs(eigenvalues[i])<5e-13)
        //     number_of_zero_eigenvalues++;
        // if (number_of_zero_eigenvalues>1)
        //   disconnected_subdomains[p] = true;
        // if (number_of_zero_eigenvalues>1)
        //   {
        //     std::cout << ">1 zero eigenvalue(s) found in subdomain " << p << std::endl;
        //     for (size_t i=0; i<eigenvalues.size(); ++i)
        //       if (std::abs(eigenvalues[i])<5e-13)
        //      std::cout << "  lambda_" << i<< " = " << eigenvalues[i] << std::endl;
        //   }

        // store the eigenvectors computed
        if (eigenvalue_fine_threshold_<=0)
          {
            // we take a fixed number of eigenvalues
            size_t n_zero=0;
            size_t n_positive=0;
            for (size_t k=0; k<eigenvalues.size(); k++) // look at all eigenvectors
              {
                if (fine_subdomain_basis_vectors[p].size()>=n_eigenvectors_fine_used_)
                  break;

                // if (eigenvalues[k]>=1e-11 && n_zero+n_positive>=n_eigenvectors_fine_used_)
                //   break; // we have enough

                // if (eigenvalues[k]<1e-11)
                //   n_zero = 1; // count them only once
                // else
                //   n_positive++;

                // rescale
                eigenvectors[k] *= 1.0/eigenvectors[k].two_norm();

                // multiply with partition of unity and zero on Dirichlet boundary
                for (size_t i=0; i<n; i++)
                  {
                    auto iglobal = pdddm->subdomainnumber_local_to_global[p][i];
                    eigenvectors[k][i] *= partition_of_unity[0][p][i]; // multiply with partition of unity, puts zero on subdomain boundary
                    for (size_t j=0; j<blocksize; j++)
                      eigenvectors[k][i][j] *= (*pglobal_dirichlet_mask)[iglobal][j]; // zero out dof on global Dirichlet boundary
                  }

                // now we have a basis vector in the subdomain
                fine_subdomain_basis_vectors[p].push_back(eigenvectors[k]);
              }
            // std::vector<size_t> colpermut(fine_subdomain_basis_vectors[p].size());
            // field_type th=1e-14;
            // field_type min_norm=1e100;
            // auto rank = rank_revealing_gram_schmidt(A,A,std::vector<VEC>(0),fine_subdomain_basis_vectors[p],colpermut,th,min_norm);
            // if (rank<fine_subdomain_basis_vectors[p].size())
            //   std::cout << "subdomain " << p << " has " << fine_subdomain_basis_vectors[p].size() << " basis vectors with " << rank << " linearly dependent" << std::endl;

            std::cout << p << ": " << fine_subdomain_basis_vectors[p].size()
                      << " basis vectors taken, first ev: " << eigenvalues[0]
                      << " last ev: " << eigenvalues[fine_subdomain_basis_vectors[p].size()-1]
                      << std::endl;
          }
        else
          {
            // adaptive choice of eigenvalues
            if (eigenvalues.size()<=1)
              DUNE_THROW(Dune::Exception, "you should probably compute more than 1 eigenvalue ...");
            size_t last_taken=0;
            for (size_t k=0; k<eigenvalues.size()-1; k++) // we know they are all nonnegative and ordered
              if (eigenvalues[k]<=eigenvalue_fine_threshold_)
                {
                  // rescale
                  eigenvectors[k] *= 1.0/eigenvectors[k].two_norm();

                  // multiply with partition of unity and zero on Dirichlet boundary
                  for (size_t i=0; i<n; i++)
                    {
                      auto iglobal = pdddm->subdomainnumber_local_to_global[p][i];
                      eigenvectors[k][i] *= partition_of_unity[0][p][i]; // multiply with partition of unity, puts zero on subdomain boundary
                      for (size_t j=0; j<blocksize; j++)
                        eigenvectors[k][i][j] *= (*pglobal_dirichlet_mask)[iglobal][j]; // zero out dof on global Dirichlet boundary
                    }

                  // now we have a basis vector in the subdomain
                  fine_subdomain_basis_vectors[p].push_back(eigenvectors[k]);
                  last_taken=k;
                }
            if (eigenvalues[last_taken+1]<=eigenvalue_fine_threshold_)
              {
                std::cout << "WARNING: number of eigenvalues computed is too small. First eigenvalue not taken: " << eigenvalues[last_taken+1] << std::endl;
              }
            std::cout << p << ": " << fine_subdomain_basis_vectors[p].size()
                      << " basis vectors taken, first ev: " << eigenvalues[0]
                      << " last ev: " << eigenvalues[fine_subdomain_basis_vectors[p].size()-1]
                      << std::endl;
          }

        // // A-orthogonalization of basis vectors
        // if (false)
        //   {
        //     std::vector<size_t> colpermut(fine_subdomain_basis_vectors[p].size());
        //     field_type min_norm;
        //     auto rank = rank_revealing_gram_schmidt(A,A,std::vector<CVEC>(0),fine_subdomain_basis_vectors[p],colpermut,abs_zero_orthogonalization_,min_norm);
        //     if (rank<fine_subdomain_basis_vectors[p].size())
        //       {
        //      std::cout << "ERROR: basis vectors are not A-orthogonal in subdomain " << p
        //                << " rank=" << rank << " size=" << fine_subdomain_basis_vectors[p].size()
        //                << std::endl;
        //      DUNE_THROW(Dune::Exception, "you should compute more than 1 eigenvalue ...");
        //       }
        //   }

        // and we set up the subdomain problem solvers for the apply phase
        // assemble dirichlet conditions on subdomain boundary
        for (size_t i=0; i<A.N(); i++)
          if (subdomain_boundary_dofs[0][p][i])
            {
              // all dofs in this block are on the boundary
              auto cIt = A[i].begin();
              auto cEndIt = A[i].end();
              for (; cIt!=cEndIt; ++cIt)
                {
                  (*cIt) = 0.0;
                  if (cIt.index()==i)
                    for (size_t comp=0; comp<blocksize; comp++) (*cIt)[comp][comp] = 1.0;
                }
            }
          else
            {
              // non dirchlet row; eliminate Dirichlet columns
              auto cIt = A[i].begin();
              auto cEndIt = A[i].end();
              for (; cIt!=cEndIt; ++cIt)
                if (subdomain_boundary_dofs[0][p][cIt.index()])
                  (*cIt) = 0.0;
            }

        // store pointer to UMFPACK solver object holding the decomposed matrix
        // hopefully it does not need the original matrix anymore ...
        timer2.reset();
#ifdef USE_CHOLMOD
        finesubdomainsolvers[p] = std::make_shared<CHOLMODSOLVER>();
        finesubdomainsolvers[p]->setMatrix(A);
#else
        finesubdomainsolvers[p] = std::make_shared<UMFPACKSOLVER>(A,false);
#endif
        setup_time_umfpack[0][p] += timer2.elapsed();
        setup_time_procs[0][p] += timer.elapsed();
      } // end loop over all subdomains on finest level

    // Now that we have produced all basis vectors we may produce the global system on the next coarser level.
    // We would not actually need that for the intermediate levels 1...Number_Of_Levels-2 but it is definitely
    // easier and allows us to do the iterated two grid method

    timer.reset();

    // first add the basisvector local<->global information to the DofMapper on this level
    std::vector<size_t> number_of_basis_vectors(nsubdomains[0]);
    for (size_t p=0; p<nsubdomains[0]; ++p)
      number_of_basis_vectors[p] = fine_subdomain_basis_vectors[p].size(); // we take from finest level
    pdddm->fill_basisvector_information(number_of_basis_vectors,disconnected_subdomains); // computes basis vector global<=> local

    // set up global coarse level matrix
    {
      auto n_coarse = pdddm->basisvector_global_to_local.size(); // now we know the size of the coarse problem
      size_t nz_coarse = 0; // compute number of nonzeroes
      for (size_t p=0; p<pdddm->p_subdomain_graph->graph.size(); p++)
        for (auto q : pdddm->p_subdomain_graph->graph[p])
          nz_coarse += pdddm->basisvector_local_to_global[p].size()*pdddm->basisvector_local_to_global[q].size();
      if (verbose_>0) std::cout << "coarse matrix " << n_coarse << "x" << n_coarse << " with " << nz_coarse << " nonzeroes" << std::endl;
      global_coarse_matrix.resize(Number_Of_Levels);
      global_coarse_matrix[1] = std::shared_ptr<CMAT>(new CMAT(n_coarse,n_coarse,nz_coarse,CMAT::row_wise)); // allocate matrix in rowwise creation mode
      for (auto row=global_coarse_matrix[1]->createbegin(); row!=global_coarse_matrix[1]->createend(); ++row) // fill nonzeroes
        {
          auto i = row.index(); // this is the global row index
          auto p = pdddm->basisvector_global_to_subdomainnumber[i]; // the subdomain this coarse vector comes from
          for (auto q : pdddm->p_subdomain_graph->graph[p]) // loop over neighboring subdomains
            for (size_t k=0; k<pdddm->basisvector_local_to_global[q].size(); k++) // loop over coarse vectors in neighboring subdomain
              {
                auto j = pdddm->basisvector_local_to_global[q][k]; // the global number of this coarse vector
                row.insert(j);
              }
        }
    }

    // assemble global coarse matrix
    for (size_t prow=0; prow<nsubdomains[0]; prow++) // loop over all subdomains (row index)
      {
        // we need to set up the subdomain matrix again because it is not stores
        std::shared_ptr<MAT> pA = set_up_local_matrix(*pAglobal,prow,pdddm->subdomainnumber_local_to_global,pdddm->subdomainnumber_global_to_local);

        // and fill it with entries of the global matrix
        // this is needed to get the correct entries on the boundary of the subdomain
        MAT& Asubdomain = *pA; // get a reference for ease of writing
        for (size_t i=0; i<Asubdomain.N(); i++) // loop over rows
          {
            auto iglobal = pdddm->subdomainnumber_local_to_global[prow][i];
            auto cIt = Asubdomain[i].begin();
            auto cEndIt = Asubdomain[i].end();
            for (; cIt!=cEndIt; ++cIt) // loop over columns
              {
                auto j = cIt.index();
                auto jglobal = pdddm->subdomainnumber_local_to_global[prow][j];
                (*cIt) = (*pAglobal)[iglobal][jglobal];
              }
          }

        // now do the triple matrix products
        VEC v1(Asubdomain.N()); // two vectors on the subdomain to compute scalar products
        VEC v2(Asubdomain.N());
        for (size_t krow=0; krow<pdddm->basisvector_local_to_global[prow].size(); krow++) // loop over all basis vectors in this subdomain
          {
            // row number in global coarse system
            auto icoarse = pdddm->basisvector_local_to_global[prow][krow];

            // matrix vector product
            Asubdomain.mv(fine_subdomain_basis_vectors[prow][krow],v1);

            // now do all the scalar products with the other basis vectors
            for (auto pcol : pdddm->p_subdomain_graph->graph[prow]) // loop over all neighboring subdomains (including myself)
              for (size_t kcol=0; kcol<pdddm->basisvector_local_to_global[pcol].size(); kcol++) // loop over all basis vectors in the neighboring subdomain
                {
                  // column number in coarse system
                  auto jcoarse = pdddm->basisvector_local_to_global[pcol][kcol];

                  // produce restriction of basis vector from subdomain pcol to subdomain prow
                  v2 = 0.0; // zero out
                  for (size_t i=0; i<Asubdomain.N(); i++)
                    {
                      auto iglobal = pdddm->subdomainnumber_local_to_global[prow][i];
                      auto it = pdddm->subdomainnumber_global_to_local[iglobal].find(pcol); // test if this global dof block exists in neighboring subdomain
                      if (it!=pdddm->subdomainnumber_global_to_local[iglobal].end()) // if yes, we have a common entry
                        {
                          auto index_in_other_subdomain = it->second;
                          v2[i] = fine_subdomain_basis_vectors[pcol][kcol][index_in_other_subdomain]; // restriction of coarse vector in subdomain pcol in subdomain prow
                        }
                    }

                  // now we have a matrix entry
                  (*global_coarse_matrix[1])[icoarse][jcoarse] = v2*v1;
                }
          }
      }

    if (verbose_>0) std::cout << "computed entries of coarse system on level " << 1 << std::endl;

    // at this point we can compute the Galerkin projection of the volume and skeleton snippets
    projected_volume_snippet_matrices.resize(Number_Of_Levels); // to store all the matrices on the different levels
    if (pskeleton_snippet_matrices!=nullptr)
      projected_skeleton_snippet_matrices.resize(Number_Of_Levels);
    if (Number_Of_Levels>2) // only then we will need projected snippet matrices
      {
        // first, project the volume snippets
        if (verbose_>0) std::cout << "projection of volume snippets to coarse basis vectors" << std::endl;
        projected_volume_snippet_matrices[0].resize(ddhierarchy[0]->number_to_volumesnippet.size()); // resize to number of snippets
        for (size_t vsn = 0; vsn<ddhierarchy[0]->number_to_volumesnippet.size(); vsn++) // loop over all volume snippets
          {
            auto& vs = ddhierarchy[0]->number_to_volumesnippet[vsn]; // get the volume snippet
            if (verbose_>1) std::cout << "Galerkin projection of " << ddhierarchy[0]->vs_to_str(vs);
            MAT& snippetA = *(*pvolume_snippet_matrices)[vsn]; // matrix for volume snippet on current level
            if (verbose_>1) std::cout << " with matrix of size " << snippetA.N() << "x" << snippetA.M();

            // we need a local to global map for the subdomain basis vectors involved with this snippet
            // this can be computed on the fly for the volume snippet
            std::map<size_t,std::map<size_t,size_t>> subdomain_local_to_snippet_local;
            size_t nc=0;
            for (auto p : vs) // loop over all subdomains intersected by the volume snippet
              for (size_t k=0; k<dddmhierarchy[0]->n_basis_vectors_per_subdomain[p]; k++)
                subdomain_local_to_snippet_local[p][k] = nc++;
            if (verbose_>1) std::cout << " nc=" << nc;

            // allocate dense matrix
            std::shared_ptr<DenseMat> pA = std::shared_ptr<DenseMat>( new DenseMat(nc,nc,0.0) ); // make a dense matrix with zeroes
            projected_volume_snippet_matrices[0][vsn] = pA; // store the matrix; it will be used on the next coarser level
            DenseMat& projectedA = *pA; // use reference

            // compute projected entries
            for (auto& colentry1 : subdomain_local_to_snippet_local) // loop over all subdomains intersected by the volume snippet
              for (auto& colentry2 : colentry1.second)
                {
                  auto pcol = colentry1.first;
                  auto kcol = colentry2.first;
                  auto jprojected = colentry2.second;

                  // restrict basis vector [pcol][kcol] to snippet
                  VEC v1(dddmhierarchy[0]->volumesnippetnumber_local_to_global[vsn].size()); v1 = 0.0;
                  for (size_t i=0; i<dddmhierarchy[0]->volumesnippetnumber_local_to_global[vsn].size(); i++) // loop over all dofs in volume snippet
                    {
                      auto iglobal = dddmhierarchy[0]->volumesnippetnumber_local_to_global[vsn][i]; // global dof
                      auto isubdomain = dddmhierarchy[0]->subdomainnumber_global_to_local[iglobal][pcol]; // local dof in subdomain pcol
                      v1[i] = fine_subdomain_basis_vectors[pcol][kcol][isubdomain]; // pick entry from basis vector
                    }

                  // multiply snippet matrix with column basis vector
                  VEC Av1(dddmhierarchy[0]->volumesnippetnumber_local_to_global[vsn].size());
                  snippetA.mv(v1,Av1);

                  // now multiply with the rows
                  for (auto& rowentry1 : subdomain_local_to_snippet_local) // loop over all subdomains intersected by the volume snippet
                    for (auto& rowentry2 : rowentry1.second)
                      {
                        auto prow = rowentry1.first;
                        auto krow = rowentry2.first;
                        auto iprojected = rowentry2.second;

                        // restrict basis vector [prow][krow] to snippet
                        VEC v2(dddmhierarchy[0]->volumesnippetnumber_local_to_global[vsn].size()); v2 = 0.0;
                        for (size_t i=0; i<dddmhierarchy[0]->volumesnippetnumber_local_to_global[vsn].size(); i++) // loop over all dofs in volume snippet
                          {
                            auto iglobal = dddmhierarchy[0]->volumesnippetnumber_local_to_global[vsn][i]; // global dof
                            auto isubdomain = dddmhierarchy[0]->subdomainnumber_global_to_local[iglobal][prow]; // local dof in subdomain prow
                            v2[i] = fine_subdomain_basis_vectors[prow][krow][isubdomain]; // pick entry from basis vector
                          }

                        // now we have one entry
                        projectedA[iprojected][jprojected] = v2 * Av1;
                      }
                }
            if (verbose_>1) std::cout << " ... done" << std::endl;
          } // loop over volume snippets

        // next are the skeleton snippets; it works basically the same ...
        if (pskeleton_snippet_matrices!=nullptr)
          {
            if (verbose_>0) std::cout << "projection of skeleton snippets to coarse basis vectors" << std::endl;
            projected_skeleton_snippet_matrices[0].resize(ddhierarchy[0]->number_to_skeletonsnippet.size()); // resize to number of snippets
            for (size_t ssn = 0; ssn<ddhierarchy[0]->number_to_skeletonsnippet.size(); ssn++) // loop over all skeleton snippets
              {
                auto& ss = ddhierarchy[0]->number_to_skeletonsnippet[ssn]; // get the skeleton snippet
                auto it = ss.begin(); // points to first v snippet number
                size_t left_volumesnippetnumber = *it;
                ++it;
                size_t right_volumesnippetnumber = *it; // so now we have the numbers of the two volume snippets to left and right
                auto& left_volumesnippet = ddhierarchy[0]->number_to_volumesnippet[left_volumesnippetnumber];
                auto& right_volumesnippet = ddhierarchy[0]->number_to_volumesnippet[right_volumesnippetnumber];

                if (verbose_>1) std::cout << "Galerkin projection of skeleton snippet " << ddhierarchy[0]->vs_to_str(ss) << " : " << ddhierarchy[0]->vs_to_str(left_volumesnippet) << "|" << ddhierarchy[0]->vs_to_str(right_volumesnippet);
                MAT& snippetA = *(*pskeleton_snippet_matrices)[ssn]; // matrix for skeleton snippet on current level
                if (verbose_>1) std::cout << " with matrix of size " << snippetA.N() << "x" << snippetA.M();

                // we need a local to global map for the subdomain basis vectors involved with this snippet
                // this can be computed on the fly for the skeleton snippet
                // NOTE: only dofs of subdomains which are BOTH left and right are involved!
                // WHY? say subdomain p is only to one side, then such a function zero (PU!) on one side and non existent on the other side
                std::map<size_t,std::map<size_t,size_t>> subdomain_local_to_snippet_local;
                size_t nc=0;
                for (auto p : left_volumesnippet) // loop over all subdomains intersected by the left volume snippet
                  if (right_volumesnippet.count(p)>0) // p is to left and right of skeleton snippet
                    for (size_t k=0; k<dddmhierarchy[0]->n_basis_vectors_per_subdomain[p]; k++)
                      subdomain_local_to_snippet_local[p][k] = nc++;
                if (verbose_>1) std::cout << " nc=" << nc;

                // allocate dense matrix
                std::shared_ptr<DenseMat> pA = std::shared_ptr<DenseMat>( new DenseMat(nc,nc,0.0) ); // make a dense matrix with zeroes
                projected_skeleton_snippet_matrices[0][ssn] = pA; // store the matrix; it will be used on the next coarser level
                DenseMat& projectedA = *pA; // use reference

                // compute projected entries
                for (auto& colentry1 : subdomain_local_to_snippet_local) // loop over all subdomains intersected by the snippet
                  for (auto& colentry2 : colentry1.second)
                    {
                      auto pcol = colentry1.first;
                      auto kcol = colentry2.first;
                      auto jprojected = colentry2.second;

                      // restrict basis vector [pcol][kcol] to snippet
                      VEC v1(dddmhierarchy[0]->skeletonsnippetnumber_local_to_global[ssn].size());
                      v1 = 0.0;
                      for (size_t i=0; i<dddmhierarchy[0]->skeletonsnippetnumber_local_to_global[ssn].size(); i++) // loop over all dofs in snippet
                        {
                          auto iglobal = dddmhierarchy[0]->skeletonsnippetnumber_local_to_global[ssn][i]; // global dof
                          auto it = dddmhierarchy[0]->subdomainnumber_global_to_local[iglobal].find(pcol);
                          if (it!=dddmhierarchy[0]->subdomainnumber_global_to_local[iglobal].end())
                            v1[i] = fine_subdomain_basis_vectors[pcol][kcol][it->second]; // pick entry from basis vector
                        }

                      // multiply snippet matrix with column basis vector
                      VEC Av1(dddmhierarchy[0]->skeletonsnippetnumber_local_to_global[ssn].size());
                      snippetA.mv(v1,Av1);

                      // now multiply with the rows
                      for (auto& rowentry1 : subdomain_local_to_snippet_local) // loop over all subdomains intersected by the volume snippet
                        for (auto& rowentry2 : rowentry1.second)
                          {
                            auto prow = rowentry1.first;
                            auto krow = rowentry2.first;
                            auto iprojected = rowentry2.second;

                            // restrict basis vector [prow][krow] to snippet
                            VEC v2(dddmhierarchy[0]->skeletonsnippetnumber_local_to_global[ssn].size());
                            v2 = 0.0;
                            for (size_t i=0; i<dddmhierarchy[0]->skeletonsnippetnumber_local_to_global[ssn].size(); i++) // loop over all dofs in snippet
                              {
                                auto iglobal = dddmhierarchy[0]->skeletonsnippetnumber_local_to_global[ssn][i]; // global dof
                                auto it = dddmhierarchy[0]->subdomainnumber_global_to_local[iglobal].find(prow);
                                if (it!=dddmhierarchy[0]->subdomainnumber_global_to_local[iglobal].end())
                                  v2[i] = fine_subdomain_basis_vectors[prow][krow][it->second]; // pick entry from basis vector
                              }

                            // now we have one entry
                            projectedA[iprojected][jprojected] = v2 * Av1;
                          }
                    }
                if (verbose_>1) std::cout << " ... done" << std::endl;
              } // loop over skeleton snippets
          } // if there are skeleton snippets
      }// if (Number_Of_Levels>2)

    setup_time_parallelizable[0] += timer.elapsed();

    //====================
    // Now we can proceed to the coarse levels; this is only needed when there are more than two levels!
    //====================

    // resize to hold level wise information
    volume_snippet_matrices.resize(Number_Of_Levels);
    if (pskeleton_snippet_matrices!=nullptr)
      skeleton_snippet_matrices.resize(Number_Of_Levels);
    coarse_subdomain_basis_vectors.resize(Number_Of_Levels);
    coarse_kerAB_vectors.resize(Number_Of_Levels);
    coarsesubdomainsolvers.resize(Number_Of_Levels);
    globalcoarsesolvers.resize(Number_Of_Levels);

    // loop over coarse levels
    for (size_t targetlevel=1; targetlevel<Number_Of_Levels; targetlevel++)
      {
        setup_time_procs[targetlevel].resize(nsubdomains[targetlevel],0.0);
        setup_time_setup_gevp[targetlevel].resize(nsubdomains[targetlevel],0.0);
        setup_time_QR[targetlevel].resize(nsubdomains[targetlevel],0.0);
        setup_time_projection[targetlevel].resize(nsubdomains[targetlevel],0.0);
        setup_time_gevp[targetlevel].resize(nsubdomains[targetlevel],0.0);
        setup_time_umfpack[targetlevel].resize(nsubdomains[targetlevel],0.0);

        std::cout << "XXX " << targetlevel << " " << setup_time_parallelizable[targetlevel] << std::endl;

        // now we entered a coarse level
        // on the finer level we already have produced the basis vectors in each subdomain
        if (verbose_>0) std::cout << "starting coarse level " << targetlevel << std::endl;

        // cluster nodes of subdomain graph on the finer level
        timer.reset();
        if (verbose_>0) std::cout << "coarse subdomain graph " << std::endl;
        auto coarsening = ddhierarchy[targetlevel-1]->p_subdomain_graph->coarsen_graph(comm,nsubdomains[targetlevel],coarse_overlap_,dddmhierarchy[targetlevel-1]->disconnected_subdomains,merge_disconnected_);
        setup_time_sequential += timer.elapsed();

        timer.reset();
        // from the coarsening we can deduce the domain decomposition on the target level
        if (verbose_>0) std::cout << "create domain decomposition information " << std::endl;
        ddhierarchy[targetlevel] = std::shared_ptr<DomainDecomposition>(new DomainDecomposition(ddhierarchy[targetlevel-1],coarsening));
        ddhierarchy[targetlevel]->print_info(verbose_);
        setup_time_parallelizable[targetlevel] += timer.elapsed();
        std::cout << "XXX after DomainDecomposition " << targetlevel << " " << setup_time_parallelizable[targetlevel] << std::endl;

        timer.reset();
        // now build dof information on this level
        if (verbose_>0) std::cout << "create dof information " << std::endl;
        dddmhierarchy[targetlevel] = std::shared_ptr<DomainDecompositionDOFMapper>(new DomainDecompositionDOFMapper(ddhierarchy[targetlevel],dddmhierarchy[targetlevel-1]) );
        dddmhierarchy[targetlevel]->print_info(verbose_);
        setup_time_parallelizable[targetlevel] += timer.elapsed();
        std::cout << "XXX after DomainDecompositionDOFMapper " << targetlevel << " " << setup_time_parallelizable[targetlevel] << std::endl;

        // if this is the coarsest level we are done; otherwise we need to compute basis vectors for the next coarser level ...
        if (targetlevel==Number_Of_Levels-1)
          {
            // factorize coarsest level system (wich already exists)
            std::cout << "factorizing coarse system on level " << targetlevel << std::endl;
            timer.reset();
            timer2.reset();
#ifdef USE_CHOLMOD
            globalcoarsesolvers[targetlevel] = std::make_shared<CCHOLMODSOLVER>();
            globalcoarsesolvers[targetlevel]->setMatrix(*global_coarse_matrix[targetlevel]);
#else
            globalcoarsesolvers[targetlevel] = std::make_shared<CUMFPACKSOLVER>(*global_coarse_matrix[targetlevel],false);
#endif
            setup_time_umfpack[targetlevel][0] += timer2.elapsed();
            setup_time_procs[targetlevel][0] += timer.elapsed();
            std::cout << "factorized coarse system on level " << targetlevel << std::endl;

            break; // exit loop over setting up coarse levels
          }

        // now here we are NOT on the coarsest level
        // extract references for ease of writing
        auto& dd = *ddhierarchy[targetlevel];
        auto& dddm = *dddmhierarchy[targetlevel];
        auto& finedd = *ddhierarchy[targetlevel-1];
        auto& finedddm = *dddmhierarchy[targetlevel-1];

        timer.reset();
        //=======================
        // Step 1: assemble volume snippet matrices on this level from projected fine grid snippet matrices
        volume_snippet_matrices[targetlevel].resize(dd.number_to_volumesnippet.size()); // space for the pointers
        if (verbose_>0) std::cout << "assemble volume snippet matrices on level " << targetlevel << std::endl;
        for (size_t cvsn=0; cvsn<dd.volumesnippet_coarse_to_fine.size(); cvsn++) // loop over all volume snippets on this level
          {
            // this is our coarse volume snippet
            auto& cvs = dd.number_to_volumesnippet[cvsn];
            if (verbose_>0) std::cout << "coarse volume snippet " << dd.vs_to_str(cvs);

            // count nonzeroes in snippet matrix
            std::map<size_t,std::set<size_t>> coupling; // coupling of the subdomains
            for (auto fvsn : dd.volumesnippet_coarse_to_fine[cvsn]) // loop over volume snippets on fine level for this volume snippet
              {
                auto& fvs = finedd.number_to_volumesnippet[fvsn];
                // all subdomains of this snippet are fully coupled
                for (auto p : fvs)
                  for (auto q : fvs)
                    coupling[p].insert(q);
              }
            if (pskeleton_snippet_matrices!=nullptr)
              for (auto fssn : dd.volumesnippet_coarse_to_fine_skeletonsnippetnumber[cvsn]) // loop over skeleton snippets on fine level contained in volume snippet
                {
                  // get snippets
                  auto& fss = finedd.number_to_skeletonsnippet[fssn];
                  auto it = fss.begin(); // points to first v snippet number
                  size_t left_volumesnippetnumber = *it;
                  ++it;
                  size_t right_volumesnippetnumber = *it; // so now we have the numbers of the two volume snippets to left and right
                  auto& left_volumesnippet = finedd.number_to_volumesnippet[left_volumesnippetnumber];
                  auto& right_volumesnippet = finedd.number_to_volumesnippet[right_volumesnippetnumber];

                  // add inter snippet couplings
                  auto s = dd.set_intersection(left_volumesnippet,right_volumesnippet);
                  for (auto p : s)
                    for (auto q : s)
                      coupling[p].insert(q);
                }

            size_t nz = 0; // compute nonzeroes from coupled subdomains
            for (auto& entry : coupling)
              for (auto& q : entry.second)
                nz += finedddm.n_basis_vectors_per_subdomain[entry.first] * finedddm.n_basis_vectors_per_subdomain[q];
            size_t n =  dddm.volumesnippetnumber_local_to_global[cvsn].size(); // number of degrees of freedom
            if (verbose_>0) std::cout << " ndof=" << n << " nz=" << nz << std::endl;

            // set up sparse matrix for volume snippet
            volume_snippet_matrices[targetlevel][cvsn] = std::shared_ptr<CMAT>(new CMAT(n,n,nz,CMAT::row_wise)); // allocate matrix in rowwise creation mode
            auto& coarsesnippetA = *volume_snippet_matrices[targetlevel][cvsn]; // reference
            for (auto row=coarsesnippetA.createbegin(); row!=coarsesnippetA.createend(); ++row) // loop over rows
              {
                auto i = row.index(); // local in volume snippet
                auto iglobal =  dddm.volumesnippetnumber_local_to_global[cvsn][i]; // global dof number on this level
                auto p = finedddm.basisvector_global_to_subdomainnumber[iglobal]; // the subdomain this coarse vector comes from
                for (auto q : coupling[p]) // loop over neighboring subdomains
                  for (size_t l=0; l<finedddm.n_basis_vectors_per_subdomain[q]; l++) // loop over basis vectors in neighboring subdomain
                    {
                      auto jglobal = finedddm.basisvector_local_to_global[q][l]; // the global number of this coarse vector
                      auto j = dddm.volumesnippetnumber_global_to_local[jglobal][cvsn]; //
                      //std::cout << "p=" << p << " q=" << q << " l=" << l << " jg=" << jglobal << " j=" << j << std::endl;
                      row.insert(j);
                    }
              }
            coarsesnippetA = 0.0; // clear matrix
            //Dune::printmatrix(std::cout, coarsesnippetA, "coarsesnippetA", "");

            // assemble from fine volume snippets contained in coarse volume snippet
            if (verbose_>0) std::cout << "assemble volume->volume " << targetlevel << std::endl;
            for (auto fvsn : dd.volumesnippet_coarse_to_fine[cvsn]) // loop over volume snippets on fine level for this volume snippet
              {
                auto& fvs = finedd.number_to_volumesnippet[fvsn];
                if (verbose_>1) std::cout << "  fine volume snippet " << dd.vs_to_str(fvs) << std::endl;

                // get projected fine level matrix for this fine snippet
                auto& finesnippetA = *projected_volume_snippet_matrices[targetlevel-1][fvsn];

                // need locol_to_local dof mapping from fine snippet to coarse snippet
                std::vector<size_t> local_to_local(finesnippetA.N());
                size_t count=0;
                for (auto p : fvs) // loop over all subdomains intersected by the volume snippet
                  for (size_t k=0; k<finedddm.n_basis_vectors_per_subdomain[p]; k++) // loop over all basis vectors in that subdomain
                    {
                      auto iglobal = finedddm.basisvector_local_to_global[p][k]; // global dof number of basis vector [p][k] on this level
                      auto ilocal_cvs = dddm.volumesnippetnumber_global_to_local[iglobal][cvsn]; // local number of that dof in coarse snippet
                      auto ilocal_fvs = count++; // local number in fine snippet
                      local_to_local[ilocal_fvs] = ilocal_cvs;
                    }

                // now assemble fine snippet to coarse snippet
                for (auto rIt=finesnippetA.begin(); rIt!=finesnippetA.end(); ++rIt)
                  {
                    auto cIt = rIt->begin();
                    auto cEndIt = rIt->end();
                    for (; cIt!=cEndIt; ++cIt)
                      {
                        //std::cout << "(" << rIt.index() << "," << cIt.index() << ") -> (" << local_to_local[rIt.index()] << "," << local_to_local[cIt.index()] << ")" << std::endl;
                        coarsesnippetA[local_to_local[rIt.index()]][local_to_local[cIt.index()]] += (*cIt);
                      }
                  }
              }

            // assemble from fine skeleton snippets contained inside the coarse volume snippet
            if (pskeleton_snippet_matrices!=nullptr)
              {
                if (verbose_>0) std::cout << "assemble skeleton->volume " << targetlevel << std::endl;
                for (auto fssn : dd.volumesnippet_coarse_to_fine_skeletonsnippetnumber[cvsn]) // loop over skeleton snippets on fine level inside this coarse volume snippet
                  {
                    auto& fss = finedd.number_to_skeletonsnippet[fssn];
                    auto it = fss.begin(); // points to first v snippet number
                    size_t left_volumesnippetnumber = *it;
                    ++it;
                    size_t right_volumesnippetnumber = *it; // so now we have the numbers of the two volume snippets to left and right
                    auto& left_volumesnippet = finedd.number_to_volumesnippet[left_volumesnippetnumber];
                    auto& right_volumesnippet = finedd.number_to_volumesnippet[right_volumesnippetnumber];
                    auto s = dd.set_intersection(left_volumesnippet,right_volumesnippet); // subdomains that are both to left and right
                    if (verbose_>1) std::cout << "  fine skeleton snippet " << finedd.ss_to_str(fss) << " intersection " << finedd.vs_to_str(s) << std::endl;

                    // get projected fine level matrix for this fine snippet
                    auto& finesnippetA = *projected_skeleton_snippet_matrices[targetlevel-1][fssn];

                    // need local_to_local dof mapping from fine snippet to coarse snippet
                    // this is computed on the fly
                    std::vector<size_t> local_to_local(finesnippetA.N());
                    size_t count=0;
                    for (auto p : left_volumesnippet) // loop over all subdomains intersected by the left volume snippet
                      if (right_volumesnippet.count(p)>0) // if p is to left and right of snippet
                        for (size_t k=0; k<finedddm.n_basis_vectors_per_subdomain[p]; k++)
                          {
                            auto iglobal = finedddm.basisvector_local_to_global[p][k]; // global dof number of basis vector [p][k] on this level
                            auto ilocal_cvs = dddm.volumesnippetnumber_global_to_local[iglobal][cvsn]; // local number of that dof in coarse snippet
                            auto ilocal_fss = count++; // local number in fine snippet
                            local_to_local[ilocal_fss] = ilocal_cvs;
                          }

                    // now assemble fine snippet to coarse snippet
                    for (auto rIt=finesnippetA.begin(); rIt!=finesnippetA.end(); ++rIt)
                      {
                        auto cIt = rIt->begin();
                        auto cEndIt = rIt->end();
                        for (; cIt!=cEndIt; ++cIt)
                          {
                            //std::cout << "(" << rIt.index() << "," << cIt.index() << ") -> (" << local_to_local[rIt.index()] << "," << local_to_local[cIt.index()] << ")" << std::endl;
                            coarsesnippetA[local_to_local[rIt.index()]][local_to_local[cIt.index()]] += (*cIt);
                          }
                      }
                  }
              }

          } // end loop over coarse volume snippets; now the volume snippet matrices are assembled;

        //=======================
        // Step 2 : assemble coarse skeleton snippets from fine skeleton snippets
        if (pskeleton_snippet_matrices!=nullptr)
          {
            skeleton_snippet_matrices[targetlevel].resize(dd.number_to_skeletonsnippet.size()); // space for the pointers
            if (verbose_>0) std::cout << "assemble skeleton snippet matrices on level " << targetlevel << std::endl;
            for (size_t cssn=0; cssn<dd.number_to_skeletonsnippet.size(); cssn++) // loop over all skelton snippets on this level
              {
                // this is our coarse volume snippet
                auto& css = dd.number_to_skeletonsnippet[cssn];
                if (verbose_>1) std::cout << "COARSE skeleton snippet " << dd.ss_to_str(css) << std::endl;

                // count nonzeroes in snippet matrix
                std::map<size_t,std::set<size_t>> coupling; // coupling of the subdomains
                for (auto fssn : dd.skeletonsnippet_coarse_to_fine[cssn]) // loop over skeleton snippets on fine level for this skeleton snippet
                  {
                    // get snippets
                    auto& fss = finedd.number_to_skeletonsnippet[fssn];
                    auto it = fss.begin(); // points to first v snippet number
                    size_t left_volumesnippetnumber = *it;
                    ++it;
                    size_t right_volumesnippetnumber = *it; // so now we have the numbers of the two volume snippets to left and right
                    auto& left_volumesnippet = finedd.number_to_volumesnippet[left_volumesnippetnumber];
                    auto& right_volumesnippet = finedd.number_to_volumesnippet[right_volumesnippetnumber];

                    // all subdomains of this snippet are fully coupled
                    auto s = finedd.set_intersection(left_volumesnippet,right_volumesnippet);
                    if (verbose_>2) std::cout << "   fine skeleton snippet " << finedd.ss_to_str(fss) << " is: " << finedd.vs_to_str(s) << std::endl;
                    for (auto p : s)
                      for (auto q : s)
                        coupling[p].insert(q);
                  }
                size_t nz = 0; // compute nonzeroes from coupled subdomains
                for (auto& entry : coupling)
                  for (auto& q : entry.second)
                    nz += finedddm.n_basis_vectors_per_subdomain[entry.first] * finedddm.n_basis_vectors_per_subdomain[q];
                size_t n =  dddm.skeletonsnippetnumber_local_to_global[cssn].size(); // number of degrees of freedom
                if (verbose_>1) std::cout << "ndof=" << n << " nz=" << nz << std::endl;

                // set up sparse matrix for skeleton snippet
                skeleton_snippet_matrices[targetlevel][cssn] = std::shared_ptr<CMAT>(new CMAT(n,n,nz,CMAT::row_wise)); // allocate matrix in rowwise creation mode
                auto& coarsesnippetA = *skeleton_snippet_matrices[targetlevel][cssn]; // reference
                for (auto row=coarsesnippetA.createbegin(); row!=coarsesnippetA.createend(); ++row) // loop over rows
                  {
                    auto i = row.index(); // local in skeleton snippet
                    auto iglobal =  dddm.skeletonsnippetnumber_local_to_global[cssn][i]; // global dof number on this level
                    auto p = finedddm.basisvector_global_to_subdomainnumber[iglobal]; // the subdomain this coarse vector comes from
                    for (auto q : coupling[p]) // loop over neighboring subdomains
                      for (size_t l=0; l<finedddm.n_basis_vectors_per_subdomain[q]; l++) // loop over basis vectors in neighboring subdomain
                        {
                          auto jglobal = finedddm.basisvector_local_to_global[q][l]; // the global number of this coarse vector
                          auto j = dddm.skeletonsnippetnumber_global_to_local[jglobal][cssn]; //
                          //std::cout << "p=" << p << " q=" << q << " l=" << l << " jg=" << jglobal << " j=" << j << std::endl;
                          row.insert(j);
                        }
                  }
                coarsesnippetA = 0.0; // clear matrix
                //Dune::printmatrix(std::cout, coarsesnippetA, "coarsesnippetA", "");

                // assemble from fine skeleton snippets
                for (auto fssn : dd.skeletonsnippet_coarse_to_fine[cssn]) // loop over skeleton snippets on fine level for this skeleton snippet
                  {
                    // get snippets
                    auto& fss = finedd.number_to_skeletonsnippet[fssn];
                    auto it = fss.begin(); // points to first v snippet number
                    size_t left_volumesnippetnumber = *it;
                    ++it;
                    size_t right_volumesnippetnumber = *it; // so now we have the numbers of the two volume snippets to left and right
                    auto& left_volumesnippet = finedd.number_to_volumesnippet[left_volumesnippetnumber];
                    auto& right_volumesnippet = finedd.number_to_volumesnippet[right_volumesnippetnumber];
                    if (verbose_>1) std::cout << "  fine skeleton snippet " << finedd.ss_to_str(fss) << std::endl;

                    // get projected fine level matrix for this fine snippet
                    auto& finesnippetA = *projected_skeleton_snippet_matrices[targetlevel-1][fssn];

                    // need local_to_local dof mapping from fine snippet to coarse snippet
                    std::vector<size_t> local_to_local;
                    for (auto p : left_volumesnippet) // loop over all subdomains intersected by the left volume snippet
                      if (right_volumesnippet.count(p)>0)
                        for (size_t k=0; k<finedddm.n_basis_vectors_per_subdomain[p]; k++)
                          {
                            auto iglobal = finedddm.basisvector_local_to_global[p][k]; // global dof number of basis vector [p][k] on this level
                            auto it =  dddm.skeletonsnippetnumber_global_to_local[iglobal].find(cssn);
                            if (it==dddm.skeletonsnippetnumber_global_to_local[iglobal].end())
                              DUNE_THROW(Dune::Exception, "assemble coarse from fine skeleton snippets: access failed");
                            auto ilocal_css = it->second; // local number of that dof in coarse snippet
                            local_to_local.push_back(ilocal_css);
                          }
                    if (local_to_local.size()!=finesnippetA.N())
                      DUNE_THROW(Dune::Exception, "assemble coarse from fine skeleton snippets: size mismatch");

                    // now assemble fine snippet to coarse snippet
                    for (auto rIt=finesnippetA.begin(); rIt!=finesnippetA.end(); ++rIt)
                      {
                        auto cIt = rIt->begin();
                        auto cEndIt = rIt->end();
                        for (; cIt!=cEndIt; ++cIt)
                          {
                            //std::cout << "(" << rIt.index() << "," << cIt.index() << ") -> (" << local_to_local[rIt.index()] << "," << local_to_local[cIt.index()] << ")" << std::endl;
                            coarsesnippetA[local_to_local[rIt.index()]][local_to_local[cIt.index()]] += (*cIt);
                          }
                      }
                  } // loop over fine snippets
              } // assembling coarse skeleton snippets from fine skeleton snippets
          } // if we have skeleton snippets at all

        // Step 3 : assemble Neumann subdomain matrices and solve eigenproblems; needs the partition of unity
        if (verbose_>0) std::cout << "assemble subdomain problems and eigensolves on level " << targetlevel << std::endl;

        // boundary dofs are determined by inspecting rows of the global matrix on this level
        subdomain_boundary_dofs[targetlevel] = dddm.subdomain_boundary;
        subdomain_boundary_mask[targetlevel] = get_subdomain_boundary_mask<field_type>(subdomain_boundary_dofs[targetlevel]);

        // build partition of unity
        partition_of_unity[targetlevel] = get_partition_of_unity_standard<field_type>(dddm.subdomainnumber_local_to_global,dddm.subdomainnumber_global_to_local,subdomain_boundary_dofs[targetlevel]);
        setup_time_parallelizable[targetlevel] += timer.elapsed();

        // loop over subdomains one by one!
        coarse_subdomain_basis_vectors[targetlevel].resize(nsubdomains[targetlevel]); // storage for basis vectors
        coarse_kerAB_vectors[targetlevel].resize(nsubdomains[targetlevel]); // storage for basis vectors
        coarsesubdomainsolvers[targetlevel].resize(nsubdomains[targetlevel]); // storage for the subdomain solvers
        std::vector<bool> disconnected_subdomains(nsubdomains[targetlevel],false); // true when subdomain is not connected
        bool negativeEigenvalues = false; // will be true if negative EVs have been found in at least one subdomain
        bool useArpack = false;
        bool useEigen = true;
        if (coarseeigensolver_=="arpack") // evaluate user flag
          {
            useArpack = true;
            useEigen = false;
          }
        for (size_t p=0; p<nsubdomains[targetlevel]; ++p)
          {
            timer.reset();

            timer2.reset();
            // set up sparse subdomain matrix
            std::shared_ptr<CMAT> pA = set_up_local_matrix(*global_coarse_matrix[targetlevel],p,dddm.subdomainnumber_local_to_global,dddm.subdomainnumber_global_to_local);
            *pA = 0.0; // clear subdomain matrix
            if (verbose_>0) std::cout << "process subdomain " << p << " on level " << targetlevel << " n=" << pA->N() << std::endl;

            // assemble snippets to subdomain matrix
            // This results in Neumann subdomain matrix without any boundary conditions!
            if (verbose_>0) std::cout << "assembling volume->subdomain" << std::endl;
            assemble_snippets_to_subdomain_matrix(p,
                                                  dd.subdomain_to_volumesnippetnumber,
                                                  dddm.volumesnippetnumber_local_to_global,
                                                  dddm.subdomainnumber_global_to_local,
                                                  volume_snippet_matrices[targetlevel],
                                                  pA);
            if (pskeleton_snippet_matrices!=nullptr)
              {
                if (verbose_>0) std::cout << "assembling interior skeleton->subdomain" << std::endl;
                assemble_snippets_to_subdomain_matrix(p,
                                                      dd.subdomain_to_interior_skeletonsnippetnumber,
                                                      dddm.skeletonsnippetnumber_local_to_global,
                                                      dddm.subdomainnumber_global_to_local,
                                                      skeleton_snippet_matrices[targetlevel],
                                                      pA);
              }
            // there are no Dirichlet conditions on the coarse levels anymore...

            // now at this point we have the NEUMANN subdomain matrix
            // now we may set up and solve the eigenvalue problems
            CMAT& A = *pA; // introduce reference for ease of writing
            make_symmetric(A); // out of paranoia
            // is_symmetric(A,"Symmetry of B",true,1e-7);

            // make the sparse B matrix out of A
            CMAT B(A);
            if (coarseGEVPrhs_=="geneo")
              {
                B = 0.0; // clear B
                assemble_snippets_to_subdomain_matrix_geneo_volume(p,
                                                                   dd.subdomain_to_volumesnippetnumber,
                                                                   dd.number_to_volumesnippet,
                                                                   dddm.volumesnippetnumber_local_to_global,
                                                                   dddm.subdomainnumber_global_to_local,
                                                                   volume_snippet_matrices[targetlevel],
                                                                   B); // assemble without interior
                if (pskeleton_snippet_matrices!=nullptr)
                  assemble_snippets_to_subdomain_matrix_geneo_skeleton(p,
                                                                       dd.subdomain_to_interior_skeletonsnippetnumber,
                                                                       dd.number_to_volumesnippet,
                                                                       dd.number_to_skeletonsnippet,
                                                                       dddm.skeletonsnippetnumber_local_to_global,
                                                                       dddm.subdomainnumber_global_to_local,
                                                                       skeleton_snippet_matrices[targetlevel],
                                                                       B);

                // std::ostringstream sBovlp;
                // sBovlp << basename << "_Bovlp_" << 1 << "_" << p << ".txt";
                // std::cout << "writing matrix file " << sBovlp.str() << std::endl;
                // writeMatrixToMatlab(B,sBovlp.str(),16);

                for (auto row_iter = B.begin(); row_iter != B.end(); ++row_iter)
                  for (auto col_iter = row_iter->begin(); col_iter != row_iter->end(); ++col_iter)
                    *col_iter *= partition_of_unity[targetlevel][p][row_iter.index()] * partition_of_unity[targetlevel][p][col_iter.index()];
              }
            else if (coarseGEVPrhs_=="1-pu")
              {
                for (auto rIt = B.begin(); rIt != B.end(); ++rIt)
                  for (auto cIt = rIt->begin(); cIt != rIt->end(); ++cIt)
                    (*cIt) *= (1.0-partition_of_unity[targetlevel][p][rIt.index()]) * (1.0-partition_of_unity[targetlevel][p][cIt.index()]);
              }
            else
              {
                for (auto rIt = B.begin(); rIt != B.end(); ++rIt)
                  for (auto cIt = rIt->begin(); cIt != rIt->end(); ++cIt)
                    (*cIt) *= partition_of_unity[targetlevel][p][rIt.index()] * partition_of_unity[targetlevel][p][cIt.index()];
              }
            setup_time_setup_gevp[targetlevel][p] = timer2.elapsed();

            // writing matrices
            // std::string basename("neumann_problem");
            // std::ostringstream sA;
            // sA << basename << "_A_" << targetlevel << "_" << p << ".txt";
            // std::cout << "writing matrix file " << sA.str() << std::endl;
            // writeMatrixToMatlab(A,sA.str(),16);
            // std::ostringstream sB;
            // sB << basename << "_B_" << targetlevel << "_" << p << ".txt";
            // std::cout << "writing matrix file " << sB.str() << std::endl;
            // writeMatrixToMatlab(B,sB.str(),16);

            //------------------------------------------------------
            // set up A+B and project it to range(A+B)
            //------------------------------------------------------

            timer2.reset();
            // counts number of basis vectors already put into the coarse space from this subdomain
            // since this can be done at various places, we define this variable here.
            size_t n_actually_used = 0;

            // QR decomposition of A+B
            size_t n=A.N(); // #dof in our coarse problem

            if (coarseeigensolver_=="straightarpack")
              {
                if (coarseGEVPrhs_=="1-pu")
                  DUNE_THROW(Dune::Exception, "straightforward application of ARPACK should only be used with partition of unity pu");
                std::cout << "straightforward application of ARPACK ... " << std::endl;

                // some additional regularization of A to make UMFPack not fail
                if (regularization_ker_>0.0)
                  for (auto rIt = A.begin(); rIt != A.end(); ++rIt)
                    for (auto cIt = rIt->begin(); cIt != rIt->end(); ++cIt)
                      if ( rIt.index()==cIt.index() && std::abs(partition_of_unity[targetlevel][p][rIt.index()])<1e-10)
                        (*cIt) += regularization_ker_;

                ArpackMLGeneo::ArPackPlusPlus_Algorithms<CMAT,CVEC> arpack(A);

                // set some parameters. Some or all of these should be input parameters.
                double tolerance = arpack_tolerance_; // tolerance for arpack algorithm
                double shift = 0.0; // shift used in arpack algorihm

                // make seperate vectors for arpack
                CVEC vec(n); // make a template
                vec = 0.0;
                std::vector<CVEC> eigenvectors(n_eigenvectors_coarse_computed_,vec); // create vector of ISTL vectors for Arpack to compute

                // make vectors to store eigenvalues. Eigenvectors will be written directly into coarse_basis_vectors.
                std::vector<double> eigenvalues(n_eigenvectors_coarse_computed_,0.0);

                // solve GenEO eigenproblems
                if (eigenvalue_coarse_threshold_<=0)
                  arpack.computeGenSymShiftInvertMinMagnitude(B,tolerance,eigenvectors,eigenvalues,shift);
                else
                  arpack.computeGenSymShiftInvertMinMagnitudeAdaptive(B,tolerance,eigenvectors,eigenvalues,shift,eigenvalue_coarse_threshold_,n_eigenvectors_coarse_used_);
                //if (verbose_>0)
                  for (size_t i=0; i<eigenvalues.size(); ++i)
                    std::cout << "  coarse lambda_" << i<< " = " << eigenvalues[i] << std::endl;

                // analyse EV for negative eigenvalues
                for (size_t i=0; i<eigenvalues.size(); ++i)
                  if (eigenvalues[i]<-1e-10)
                    {
                      CVEC Ae(n);
                      A.mv(eigenvectors[i],Ae);
                      CVEC Be(n);
                      B.mv(eigenvectors[i],Be);
                      auto X(Ae); X.axpy(-eigenvalues[i],Be);
                      auto Y(Ae); Y.axpy(eigenvalues[i],Be);
                      std::cout << i << ": " << "||e||=" << eigenvectors[i].infinity_norm()
                                << " ||Ae||=" << Ae.infinity_norm()
                                << " ||Be||=" << Be.infinity_norm()
                                << " ||Ae-lambda*Be||=" << X.infinity_norm()
                                << " ||Ae+lambda*Be||=" << Y.infinity_norm()
                                << std::endl;
                    }

                // store eigenvectors
                for (size_t k=0; k<eigenvalues.size(); k++)
                  if ( (eigenvalue_coarse_threshold_<=0 && n_actually_used<n_eigenvectors_coarse_used_) || (eigenvalue_coarse_threshold_>0 && eigenvalues[k]<=eigenvalue_coarse_threshold_) )
                    if (eigenvalues[k]>-1e-12)
                      {
                        // generate eigenvector
                        CVEC& v = eigenvectors[k];

                        // multiply with partition of unity
                        for (size_t i=0; i<n; i++)
                          v[i] *= partition_of_unity[targetlevel][p][i]; // multiply with partition of unity, puts zero on subdomain boundary

                        // rescale
                        if (verbose_>1) std::cout << "norm of basisvector " << k << " is " << v.two_norm() << std::endl;
                        v *= 1.0/v.two_norm();

                        // now we have a basis vector in the subdomain
                        coarse_subdomain_basis_vectors[targetlevel][p].push_back(v);
                        n_actually_used++;
                      }
                std::cout << p << ": " << coarse_subdomain_basis_vectors[targetlevel][p].size() << " basis vectors taken" << std::endl;

                // check adaptive ev threshold
                if (eigenvalue_coarse_threshold_>0)
                  {
                    if (n_actually_used==eigenvalues.size())
                      std::cout << "WARNING: all eigenvectors used, cannot guarantee threshold" << std::endl;
                    else
                      if (eigenvalues[n_actually_used]<=eigenvalue_coarse_threshold_)
                        std::cout << "WARNING: number of eigenvalues not sufficient for threshold. First not taken: " << eigenvalues[n_actually_used] << std::endl;
                  }
                setup_time_gevp[targetlevel][p] = timer2.elapsed();
              }
            else
              {
                // Version with QR decomposition of A+B to construct orthogonal complement of ker(A) cap ker(B)
                // make A+B as eigen matrix
                using EM = Eigen::Matrix<field_type,Eigen::Dynamic,Eigen::Dynamic>;
                using EV = Eigen::Matrix<field_type,Eigen::Dynamic,1>;
                std::vector<CVEC> kerA_intersect_kerB_complement;
                {
                  EM AB(n,n);
                  for (size_t j=0; j<n; j++)
                    for (size_t i=0; i<n; i++)
                      AB(i,j) = 0.0;
                  for (auto rIt = A.begin(); rIt != A.end(); ++rIt)
                    for (auto cIt = rIt->begin(); cIt != rIt->end(); ++cIt)
                      AB(rIt.index(),cIt.index()) = (*cIt);
                  for (auto rIt = B.begin(); rIt != B.end(); ++rIt)
                    for (auto cIt = rIt->begin(); cIt != rIt->end(); ++cIt)
                      AB(rIt.index(),cIt.index()) += (*cIt);
                  if (regularization_ker_>0.0)
                    for (size_t i=0; i<n; i++)
                      AB(i,i) += regularization_ker_;

                  // compute QR decomposition
                  Eigen::ColPivHouseholderQR<EM> qr;
                  if (abs_zero_ker_>0.0) // else take the default threshold
                    qr.setThreshold(abs_zero_ker_*AB.template lpNorm<Eigen::Infinity>());
                  std::cout << p << ": QR decomposition of A+B of size " << n << " ... ";
                  qr.compute(AB);
                  auto rankAB = qr.rank();
                  std::cout << " done. rank=" << rankAB << " nonzero pivots=" << qr.nonzeroPivots() << std::endl;
                  auto householder = qr.householderQ();
                  householder.setLength(rankAB);
                  EM Q = householder;

                  // the first rankAB columns of Q form an orthogonal basis of orthogonal complement of ker(A) cap ker(B)
                  CVEC v(n);
                  for (size_t j=0; j<rankAB; j++)
                    {
                      for (size_t i=0; i<n; i++)
                        v[i] = Q(i,j); // copy column j
                      kerA_intersect_kerB_complement.push_back(v); // add column
                    }

                  // debugging code to visualize the kernel of A+B
                  if (false)
                    {
                      EM R = qr.matrixR(); // .topLeftCorner(qr.rank(), qr.rank()).template triangularView<Eigen::Upper>();
                      for (size_t k=rankAB; k<n; k++)
                        {
                          EV w(n);
                          for (auto i=0; i<n; i++) w(i) = 0.0;
                          w(k) = 1.0;
                          // backsolve
                          for (int i=rankAB-1; i>=0; i--)
                            {
                              field_type sum = 0.0;
                              for (auto j=i+1; j<n; j++)
                                sum += R(i,j)*w(j);
                              w(i) = -sum/R(i,i);
                            }
                          EV r(n);
                          r = qr.colsPermutation()*w; // permutation

                          // copy to Dune Vector
                          CVEC v(n);
                          for (auto i=0; i<n; i++)
                            v[i] = r(i); // was: Q(i,k);

                          coarse_kerAB_vectors[targetlevel][p].push_back(v);
                        }
                    }// end computation of ker(A) \cap ker(B) and \ker(A) without \ker(B)
                } // QR decomposition
                setup_time_QR[targetlevel][p] = timer2.elapsed();

                // project eigenvalue problem to the orthogonal complement of ker(a) cap ker(B)
                // std::cout << "project eigenvalue problem to subspace" << std::endl;

                // the size of the projected problem
                auto N = kerA_intersect_kerB_complement.size();

                if (useEigen)
                  {
                    timer2.reset();
                    // project A to subspace
                    EM projectedA(N,N);
                    for (size_t j=0; j<N; j++)
                      {
                        CVEC Av(n);
                        A.mv(kerA_intersect_kerB_complement[j],Av);
                        for (size_t i=0; i<N; i++)
                          projectedA(i,j) = kerA_intersect_kerB_complement[i]*Av;
                      }

                    // project B to subspace
                    EM projectedB(N,N);
                    for (size_t j=0; j<N; j++)
                      {
                        CVEC Bv(n);
                        B.mv(kerA_intersect_kerB_complement[j],Bv);
                        for (size_t i=0; i<N; i++)
                          projectedB(i,j) = kerA_intersect_kerB_complement[i]*Bv;
                      }
                    setup_time_projection[targetlevel][p] = timer2.elapsed();

                    // solve eigenproblem
                    std::cout << p << ": solve eigenvalue problem using eigen ... ";
                    timer2.reset();
                    Eigen::GeneralizedEigenSolver<EM> ges;
                    ges.compute(projectedA,projectedB);
                    std::cout << "done" << std::endl;

                    // and get the result
                    auto alphas = ges.alphas();
                    auto betas = ges.betas();
                    auto evs = ges.eigenvectors();
                    // analyse and sort eigenvalues
                    std::vector<std::pair<field_type,size_t>> sortedevs;
                    for (int k=0; k<alphas.size(); k++)
                      {
                        //std::cout << k << ": alpha=" << alphas(k) << " beta=" << betas(k) << std::endl;
                        if (std::abs(alphas(k).imag())>0.25*std::abs(alphas(k).real()))
                          {
                            std::cout << k << ": imaginary eigenvalue alpha=" << alphas(k) << " beta=" << betas(k) << " skipping" << std::endl;
                            //DUNE_THROW(Dune::Exception, "imaginary eigenvalue found");
                          }
                        else
                          if (betas(k)!=0.0)
                            if (alphas(k).real()/betas(k)>-1e-5 && alphas(k).real()/betas(k)<1e6)
                              sortedevs.push_back(std::pair<field_type,size_t>(alphas(k).real()/betas(k),k));
                            else
                              if (std::abs(alphas(k).real()/betas(k))<1e6 && alphas(k).real()/betas(k)<=-1e-5)
                                std::cout << "EV negative: k=" << k << " alpha=" << alphas(k) << " beta=" << betas(k) << " lambda=" << alphas(k).real()/betas(k) << std::endl;
                      }
                    if (sortedevs.size()==0)
                      {
                        std::cout << "Ups, sort container is zero ! " << alphas.size() << std::endl;
                        for (int k=0; k<alphas.size(); k++)
                          std::cout << k << ": alpha=" << alphas(k) << " beta=" << betas(k) << std::endl;
                        exit(1);
                      }
                    std::sort(sortedevs.begin(),sortedevs.end(),[&](const std::pair<field_type,size_t>& a, const std::pair<field_type,size_t>& b) {return a.first<b.first;});
                    if (sortedevs.size()>0)
                      if (verbose_>1)
                        for (int i=0; i<std::min(sortedevs.size(),n_eigenvectors_coarse_computed_); i++)
                          std::cout << "lambda " << i << " (" << sortedevs[i].second << ") computed with eigen: " << sortedevs[i].first << std::endl;

                    // store eigenvectors
                    size_t last_k=0;
                    for (size_t k=0; k<std::min(sortedevs.size(),n_eigenvectors_coarse_computed_); k++)
                      if ( (eigenvalue_coarse_threshold_<=0 && n_actually_used<n_eigenvectors_coarse_used_) || (eigenvalue_coarse_threshold_>0 && sortedevs[k].first<=eigenvalue_coarse_threshold_) )
                        if ( sortedevs[k].first>-1e-10 )
                          {
                            // generate eigenvector
                            CVEC v(A.N());
                            v = 0;
                            size_t j = sortedevs[k].second; // original number of ev before sorting
                            for (size_t i=0; i<N; i++)
                              v.axpy(evs(i,j).real(),kerA_intersect_kerB_complement[i]);

                            // check how well this is an eigenvector
                            if (verbose_>1)
                              {
                                CVEC Av(n);
                                A.mv(v,Av);
                                CVEC Bv(n);
                                B.mv(v,Bv);
                                CVEC z(Av);
                                z.axpy(-sortedevs[k].first,Bv);
                                std::cout << k << ": EIG adding basis vector |v|=" << v.two_norm() << " |Av|=" << Av.two_norm() << " |Bv|=" << Bv.two_norm() << " |Av-lambda*Bv|=" << z.two_norm() << std::endl;
                              }

                            // multiply with partition of unity
                            for (size_t i=0; i<A.N(); i++)
                              v[i] *= partition_of_unity[targetlevel][p][i]; // multiply with partition of unity, puts zero on subdomain boundary

                            // rescale
                            v *= 1.0/v.two_norm();

                            // now we have a basis vector in the subdomain
                            coarse_subdomain_basis_vectors[targetlevel][p].push_back(v);
                            n_actually_used++;
                            last_k = k;
                          }

                    std::cout << p << ": " << n_actually_used << " basis vectors taken, last ev: " << sortedevs[last_k].first << std::endl;
                    if (n_actually_used==0)
                      DUNE_THROW(Dune::Exception, "no eigenvectors taken in adaptive choice of eigenvectors");

                    // check adaptive ev threshold
                    if (eigenvalue_coarse_threshold_>0)
                      {
                        if (n_actually_used==sortedevs.size())
                          std::cout << "WARNING: all eigenvectors used, cannot guarantee threshold" << std::endl;
                        else
                          if (sortedevs[n_actually_used+1].first<=eigenvalue_coarse_threshold_)
                            std::cout << "WARNING: number of eigenvalues not sufficient for threshold. First not taken: " << sortedevs[n_actually_used+1].first << std::endl;
                      }
                    setup_time_gevp[targetlevel][p] = timer2.elapsed();
                  }

                // solve the EVP with Arpack
                if (useArpack)
                  {
                    timer2.reset();
                    std::cout << "projection of EVP to subspace ... " << std::endl;
                    // now we need to set up the projected matrices
                    // the Arpack Wrapper can only handle sparse matrices, but ours are full :-(
                    CMAT projectedA(N,N,N*N,CMAT::row_wise);
                    for (auto row=projectedA.createbegin(); row!=projectedA.createend(); ++row) // fill nonzeroes
                      for (size_t j=0; j<N; j++)
                        row.insert(j);
                    //std::cout << "A:";
                    for (size_t j=0; j<N; j++)
                      {
                        CVEC v(A.N());
                        A.mv(kerA_intersect_kerB_complement[j],v);
                        //std::cout << j << "|" << v.two_norm() << " ";
                        for (size_t i=0; i<N; i++)
                          projectedA[i][j] = kerA_intersect_kerB_complement[i]*v;
                      }
                    make_symmetric(projectedA);
                    //std::cout << std::endl;
                    //Dune::printmatrix(std::cout, projectedA, "projectedA", "");

                    CMAT projectedB(projectedA); // same structure as A
                    //std::cout << "B:";
                    for (size_t j=0; j<N; j++)
                      {
                        CVEC v(B.N());
                        B.mv(kerA_intersect_kerB_complement[j],v);
                        //std::cout << j << "|" << v.two_norm() << " ";
                        for (size_t i=0; i<N; i++)
                          projectedB[i][j] = kerA_intersect_kerB_complement[i]*v;
                      }
                    make_symmetric(projectedB);
                    std::cout << "done" << std::endl;
                    //Dune::printmatrix(std::cout, projectedB, "projectedB", "");
                    setup_time_projection[targetlevel][p] = timer2.elapsed();

                    std::cout << "solve eigenvalue problem using arpack ... ";
                    timer2.reset();
                    // now solve the projected eigenproblem
                    // Setup ArpackGenEO wrapper, solving a generalized eigenproblem with lhs A.
                    ArpackMLGeneo::ArPackPlusPlus_Algorithms<CMAT,CVEC> projected_arpack(projectedA);

                    // set some parameters. Some or all of these should be input parameters.
                    double tolerance = arpack_tolerance_; // tolerance for arpack algorithm
                    double shift = 0.001; // shift used in arpack algorihm

                    // make seperate vectors for arpack
                    CVEC vec(N); // make a template
                    vec = 0.0;
                    std::vector<CVEC> eigenvectors(n_eigenvectors_coarse_computed_,vec); // create vector of ISTL vectors for Arpack to compute

                    // make vectors to store eigenvalues. Eigenvectors will be written directly into coarse_basis_vectors.
                    std::vector<double> eigenvalues(n_eigenvectors_coarse_computed_,0.0);

                    // solve GenEO eigenproblems
                    projected_arpack.computeStdNonSymMinMagnitude(projectedB,tolerance,eigenvectors,eigenvalues,shift);
                    if (verbose_>0)
                      for (size_t i=0; i<eigenvalues.size(); ++i)
                        std::cout << "  projected lambda_" << i<< " = " << eigenvalues[i] << std::endl;

                    // analyse EV for negative eigenvalues
                    for (size_t i=0; i<eigenvalues.size(); ++i)
                      if (eigenvalues[i]<-1e-10)
                        {
                          CVEC Ae(N);
                          A.mv(eigenvectors[i],Ae);
                          CVEC Be(N);
                          B.mv(eigenvectors[i],Be);
                          auto X(Ae); X.axpy(-eigenvalues[i],Be);
                          auto Y(Ae); Y.axpy(eigenvalues[i],Be);
                          std::cout << i << ": " << "||e||=" << eigenvectors[i].infinity_norm()
                                    << " ||Ae||=" << Ae.infinity_norm()
                                    << " ||Be||=" << Be.infinity_norm()
                                    << " ||Ae-lambda*Be||=" << X.infinity_norm()
                                    << " ||Ae+lambda*Be||=" << Y.infinity_norm()
                                    << std::endl;
                        }

                    // store eigenvectors
                    for (size_t k=0; k<eigenvalues.size(); k++)
                      if ( (eigenvalue_coarse_threshold_<=0 && n_actually_used<n_eigenvectors_coarse_used_) || (eigenvalue_coarse_threshold_>0 && eigenvalues[k]<=eigenvalue_coarse_threshold_) )
                        if (eigenvalues[k]>-1e-12)
                          {
                            // generate eigenvector
                            CVEC v(A.N());
                            v = 0;
                            for (size_t i=0; i<N; i++)
                              v.axpy(eigenvectors[k][i],kerA_intersect_kerB_complement[i]);

                            // check how well this is an eigenvector
                            if (verbose_>1)
                              {
                                CVEC Av(n);
                                A.mv(v,Av);
                                CVEC Bv(n);
                                B.mv(v,Bv);
                                CVEC z(Av);
                                z.axpy(-eigenvalues[k],Bv);
                                std::cout << k << ": ARP adding basis vector |v|=" << v.two_norm() << " |Av|=" << Av.two_norm() << " |Bv|=" << Bv.two_norm() << " |Av-lambda*Bv|=" << z.two_norm() << std::endl;
                              }

                            // multiply with partition of unity
                            for (size_t i=0; i<A.N(); i++)
                              v[i] *= partition_of_unity[targetlevel][p][i]; // multiply with partition of unity, puts zero on subdomain boundary

                            // rescale
                            if (verbose_>1) std::cout << "norm of basisvector " << k << " is " << v.two_norm() << std::endl;
                            v *= 1.0/v.two_norm();

                            // now we have a basis vector in the subdomain
                            coarse_subdomain_basis_vectors[targetlevel][p].push_back(v);
                            n_actually_used++;
                          }
                    std::cout << p << ": " << coarse_subdomain_basis_vectors[targetlevel][p].size() << " basis vectors taken" << std::endl;

                    // check adaptive ev threshold
                    if (eigenvalue_coarse_threshold_>0)
                      {
                        if (n_actually_used==eigenvalues.size())
                          std::cout << "WARNING: all eigenvectors used, cannot guarantee threshold" << std::endl;
                        else
                          if (eigenvalues[n_actually_used+1]<=eigenvalue_coarse_threshold_)
                            std::cout << "WARNING: number of eigenvalues not sufficient for threshold. First not taken: " << eigenvalues[n_actually_used+1] << std::endl;
                      }
                    setup_time_gevp[targetlevel][p] = timer2.elapsed();
                  }
              } // solve projected GEVP

            //------------------------------------------------------
            // now make the subdomain matrix out of A and factorize it
            //------------------------------------------------------

            // assemble dirichlet conditions on subdomain boundary in A
            for (size_t i=0; i<A.N(); i++)
              if (subdomain_boundary_dofs[targetlevel][p][i])
                {
                  // boundary row
                  auto cIt = A[i].begin();
                  auto cEndIt = A[i].end();
                  for (; cIt!=cEndIt; ++cIt)
                    (*cIt) = (cIt.index()==i) ? 1.0 : 0.0;
                }
              else
                {
                  // non boundary row
                  auto cIt = A[i].begin();
                  auto cEndIt = A[i].end();
                  for (; cIt!=cEndIt; ++cIt)
                    if (subdomain_boundary_dofs[targetlevel][p][cIt.index()])
                      (*cIt) = 0.0;
                }

            // // compare with C
            // is_equal(A,C,"A=C",true,1e-12);

            // factorize subdomain system with UMFPack
            timer2.reset();
#ifdef USE_CHOLMOD
            coarsesubdomainsolvers[targetlevel][p] = std::make_shared<CCHOLMODSOLVER>();
            coarsesubdomainsolvers[targetlevel][p]->setMatrix(A);
#else
            coarsesubdomainsolvers[targetlevel][p] = std::make_shared<CUMFPACKSOLVER>(A,false);
#endif
            setup_time_umfpack[targetlevel][p] += timer2.elapsed();

            setup_time_procs[targetlevel][p] += timer.elapsed();
          }

        // if (negativeEigenvalues)
        //   DUNE_THROW(Dune::Exception, "negative eigenvalue found");

        timer.reset();

        // as result there are now basis vectors on this level
        std::vector<size_t> number_of_basis_vectors(nsubdomains[targetlevel]);
        for (size_t p=0; p<nsubdomains[targetlevel]; p++)
          number_of_basis_vectors[p] = coarse_subdomain_basis_vectors[targetlevel][p].size();
        dddm.fill_basisvector_information(number_of_basis_vectors,disconnected_subdomains); // computes basis vector global<=> local

        // and we need to construct the global coarse matrix ...
        // set up global coarse level matrix
        {
          auto n_coarse = dddm.basisvector_global_to_local.size(); // now we know the size of the coarse problem
          size_t nz_coarse = 0; // compute number of nonzeroes
          for (size_t p=0; p<dddm.p_subdomain_graph->graph.size(); p++)
            for (auto q : dddm.p_subdomain_graph->graph[p])
              nz_coarse += dddm.basisvector_local_to_global[p].size()*dddm.basisvector_local_to_global[q].size();
          if (verbose_>0) std::cout << "level " << targetlevel+1 << " coarse matrix " << n_coarse << "x" << n_coarse << " with " << nz_coarse << " nonzeroes" << std::endl;
          global_coarse_matrix[targetlevel+1] = std::shared_ptr<CMAT>(new CMAT(n_coarse,n_coarse,nz_coarse,CMAT::row_wise)); // allocate matrix in rowwise creation mode
          CMAT& matrix = *global_coarse_matrix[targetlevel+1];
          for (auto row=matrix.createbegin(); row!=matrix.createend(); ++row) // fill nonzeroes
            {
              auto i = row.index(); // this is the global row index
              auto p = dddm.basisvector_global_to_subdomainnumber[i]; // the subdomain this coarse vector comes from
              for (auto q : dddm.p_subdomain_graph->graph[p]) // loop over neighboring subdomains
                for (size_t k=0; k<dddm.basisvector_local_to_global[q].size(); k++) // loop over coarse vectors in neighboring subdomain
                  {
                    auto j = dddm.basisvector_local_to_global[q][k]; // the global number of this coarse vector
                    row.insert(j);
                  }
            }
        }

        // assemble global coarse matrix
        for (size_t p_col=0; p_col<nsubdomains[targetlevel]; p_col++) // loop over all subdomains (row index) on this level; each subdomain gives rise to a block row on coarse level
          {
            // we need to set up the subdomain matrix again
            std::shared_ptr<CMAT> pA = set_up_local_matrix(*global_coarse_matrix[targetlevel],p_col,dddm.subdomainnumber_local_to_global,dddm.subdomainnumber_global_to_local);

            // ... by filling it with entries of the global matrix
            // this is needed to get the correct entries on the boundary of the subdomain (there, the Neumann matrix would be wrong)
            CMAT& Asubdomain = *pA; // get a reference for ease of writing
            for (size_t i=0; i<Asubdomain.N(); i++) // loop over rows
              {
                auto iglobal = dddm.subdomainnumber_local_to_global[p_col][i];
                auto cIt = Asubdomain[i].begin();
                auto cEndIt = Asubdomain[i].end();
                for (; cIt!=cEndIt; ++cIt) // loop over columns
                  {
                    auto j = cIt.index();
                    auto jglobal = dddm.subdomainnumber_local_to_global[p_col][j];
                    (*cIt) = (*global_coarse_matrix[targetlevel])[iglobal][jglobal];
                  }
              }

            // now do the triple matrix products
            CVEC v1(Asubdomain.N()); // two vectors on the subdomain to compute scalar products
            CVEC v2(Asubdomain.N());
            for (size_t k_col=0; k_col<dddm.basisvector_local_to_global[p_col].size(); k_col++) // loop over all basis vectors in this subdomain
              {
                // column number in global coarse system
                auto j = dddm.basisvector_local_to_global[p_col][k_col];

                // matrix vector product A*r_j
                Asubdomain.mv(coarse_subdomain_basis_vectors[targetlevel][p_col][k_col],v1);

                // now do all the scalar products with the other (row) basis vectors
                for (auto p_row : dddm.p_subdomain_graph->graph[p_col]) // loop over all neighboring subdomains (including myself)
                  for (size_t k_row=0; k_row<dddm.basisvector_local_to_global[p_row].size(); k_row++) // loop over all basis vectors in the neighboring subdomain
                    {
                      // row number in coarse system
                      auto i = dddm.basisvector_local_to_global[p_row][k_row];

                      // produce restriction of basis vector from subdomain p_row to subdomain p_col
                      v2 = 0.0; // zero out
                      for (size_t i=0; i<Asubdomain.N(); i++)
                        {
                          auto iglobal = dddm.subdomainnumber_local_to_global[p_col][i];
                          auto it = dddm.subdomainnumber_global_to_local[iglobal].find(p_row); // test if this global dof block exists in neighboring subdomain
                          if (it!=dddm.subdomainnumber_global_to_local[iglobal].end()) // if yes, we have a common entry
                            {
                              auto index_in_other_subdomain = it->second;
                              v2[i] = coarse_subdomain_basis_vectors[targetlevel][p_row][k_row][index_in_other_subdomain]; // restriction of coarse vector in subdomain p_row in subdomain p_col
                            }
                        }

                      // now we have a matrix entry
                      (*global_coarse_matrix[targetlevel+1])[i][j] = v2*v1;
                    }
              }
          }
        setup_time_parallelizable[targetlevel] += timer.elapsed();
        if (verbose_>0) std::cout << "computed entries of coarse system on level " << targetlevel+1 << std::endl;

        timer.reset();

        // at this point we can compute the Galerkin projection of the volume and skeleton snippets to the coarse basis vectors
        if (targetlevel<Number_Of_Levels-2) // if next coarse level is not the coarsest
          {
            // first, project the volume snippets
            if (verbose_>0) std::cout << "projection of volume snippets to coarse basis vectors" << std::endl;
            projected_volume_snippet_matrices[targetlevel].resize(dd.number_to_volumesnippet.size()); // resize to number of snippets
            for (size_t vsn = 0; vsn<dd.number_to_volumesnippet.size(); vsn++) // loop over all volume snippets on this level
              {
                auto& vs = dd.number_to_volumesnippet[vsn]; // get the volume snippet
                if (verbose_>1) std::cout << "Galerkin projection of " << dd.vs_to_str(vs);
                CMAT& snippetA = *volume_snippet_matrices[targetlevel][vsn]; // matrix for volume snippet on current level
                if (verbose_>1) std::cout << " with matrix of size " << snippetA.N() << "x" << snippetA.M();

                // we need a local to global map for the subdomain basis vectors involved with this snippet
                // this can be computed on the fly for the volume snippet
                std::map<size_t,std::map<size_t,size_t>> subdomain_local_to_snippet_local;
                size_t nc=0;
                for (auto p : vs) // loop over all subdomains intersected by the volume snippet
                  for (size_t k=0; k<dddm.n_basis_vectors_per_subdomain[p]; k++)
                    subdomain_local_to_snippet_local[p][k] = nc++;
                if (verbose_>1) std::cout << " nc=" << nc;

                // allocate dense matrix
                std::shared_ptr<DenseMat> pA = std::shared_ptr<DenseMat>( new DenseMat(nc,nc,0.0) ); // make a dense matrix with zeroes
                projected_volume_snippet_matrices[targetlevel][vsn] = pA; // store the matrix; it will be used on the next coarser level
                DenseMat& projectedA = *pA; // use reference

                // compute projected entries
                for (auto& colentry1 : subdomain_local_to_snippet_local) // loop over all subdomains intersected by the volume snippet
                  for (auto& colentry2 : colentry1.second)
                    {
                      auto pcol = colentry1.first;
                      auto kcol = colentry2.first;
                      auto jprojected = colentry2.second;

                      // restrict basis vector [pcol][kcol] to snippet
                      CVEC v1(dddm.volumesnippetnumber_local_to_global[vsn].size()); v1 = 0.0;
                      for (size_t i=0; i<dddm.volumesnippetnumber_local_to_global[vsn].size(); i++) // loop over all dofs in volume snippet
                        {
                          auto iglobal = dddm.volumesnippetnumber_local_to_global[vsn][i]; // global dof
                          auto isubdomain = dddm.subdomainnumber_global_to_local[iglobal][pcol]; // local dof in subdomain pcol
                          v1[i] = coarse_subdomain_basis_vectors[targetlevel][pcol][kcol][isubdomain]; // pick entry from basis vector
                        }

                      // multiply snippet matrix with column basis vector
                      CVEC Av1(dddm.volumesnippetnumber_local_to_global[vsn].size());
                      snippetA.mv(v1,Av1);

                      // now multiply with the rows
                      for (auto& rowentry1 : subdomain_local_to_snippet_local) // loop over all subdomains intersected by the volume snippet
                        for (auto& rowentry2 : rowentry1.second)
                          {
                            auto prow = rowentry1.first;
                            auto krow = rowentry2.first;
                            auto iprojected = rowentry2.second;

                            // restrict basis vector [prow][krow] to snippet
                            CVEC v2(dddm.volumesnippetnumber_local_to_global[vsn].size()); v2 = 0.0;
                            for (size_t i=0; i<dddm.volumesnippetnumber_local_to_global[vsn].size(); i++) // loop over all dofs in volume snippet
                              {
                                auto iglobal = dddm.volumesnippetnumber_local_to_global[vsn][i]; // global dof
                                auto isubdomain = dddm.subdomainnumber_global_to_local[iglobal][prow]; // local dof in subdomain prow
                                v2[i] = coarse_subdomain_basis_vectors[targetlevel][prow][krow][isubdomain]; // pick entry from basis vector
                              }

                            // now we have one entry
                            projectedA[iprojected][jprojected] = v2 * Av1;
                          }
                    }
                if (verbose_>1) std::cout << " ... done" << std::endl;
              } // loop over volume snippets

                // next are the skeleton snippets; it works basically the same ...
            if (pskeleton_snippet_matrices!=nullptr)
              {
                if (verbose_>0) std::cout << "projection of skeleton snippets to coarse basis vectors" << std::endl;
                projected_skeleton_snippet_matrices[targetlevel].resize(dd.number_to_skeletonsnippet.size()); // resize to number of snippets
                for (size_t ssn = 0; ssn<dd.number_to_skeletonsnippet.size(); ssn++) // loop over all skeleton snippets
                  {
                    auto& ss = dd.number_to_skeletonsnippet[ssn]; // get the skeleton snippet
                    auto it = ss.begin(); // points to first v snippet number
                    size_t left_volumesnippetnumber = *it;
                    ++it;
                    size_t right_volumesnippetnumber = *it; // so now we have the numbers of the two volume snippets to left and right
                    auto& left_volumesnippet = dd.number_to_volumesnippet[left_volumesnippetnumber];
                    auto& right_volumesnippet = dd.number_to_volumesnippet[right_volumesnippetnumber];

                    if (verbose_>1) std::cout << "Galerkin projection of skeleton snippet " << dd.vs_to_str(ss) << " : " << dd.vs_to_str(left_volumesnippet) << "|" << dd.vs_to_str(right_volumesnippet);
                    CMAT& snippetA = *skeleton_snippet_matrices[targetlevel][ssn]; // matrix for skeleton snippet on current level
                    if (verbose_>1) std::cout << " with matrix of size " << snippetA.N() << "x" << snippetA.M();

                    // we need a local to global map for the subdomain basis vectors involved with this snippet
                    // this can be computed on the fly for the skeleton snippet
                    // NOTE: only dofs of subdomains which are BOTH left and right are involved!
                    // WHY? say subdomain p is only to one side, then such a function zero (PU!) on one side and non existent on the other side
                    std::map<size_t,std::map<size_t,size_t>> subdomain_local_to_snippet_local;
                    size_t nc=0;
                    for (auto p : left_volumesnippet) // loop over all subdomains intersected by the left volume snippet
                      if (right_volumesnippet.count(p)>0) // p is to left and right of skeleton snippet
                        for (size_t k=0; k<dddm.n_basis_vectors_per_subdomain[p]; k++)
                          subdomain_local_to_snippet_local[p][k] = nc++;
                    if (verbose_>1) std::cout << " nc=" << nc;

                    // allocate dense matrix
                    std::shared_ptr<DenseMat> pA = std::shared_ptr<DenseMat>( new DenseMat(nc,nc,0.0) ); // make a dense matrix with zeroes
                    projected_skeleton_snippet_matrices[targetlevel][ssn] = pA; // store the matrix; it will be used on the next coarser level
                    DenseMat& projectedA = *pA; // use reference

                    // compute projected entries
                    for (auto& colentry1 : subdomain_local_to_snippet_local) // loop over all subdomains intersected by the snippet
                      for (auto& colentry2 : colentry1.second)
                        {
                          auto pcol = colentry1.first;
                          auto kcol = colentry2.first;
                          auto jprojected = colentry2.second;

                          // restrict basis vector [pcol][kcol] to snippet
                          CVEC v1(dddm.skeletonsnippetnumber_local_to_global[ssn].size());
                          v1 = 0.0;
                          for (size_t i=0; i<dddm.skeletonsnippetnumber_local_to_global[ssn].size(); i++) // loop over all dofs in snippet
                            {
                              auto iglobal = dddm.skeletonsnippetnumber_local_to_global[ssn][i]; // global dof
                              auto it = dddm.subdomainnumber_global_to_local[iglobal].find(pcol);
                              if (it!=dddm.subdomainnumber_global_to_local[iglobal].end())
                                v1[i] = coarse_subdomain_basis_vectors[targetlevel][pcol][kcol][it->second]; // pick entry from basis vector
                            }

                          // multiply snippet matrix with column basis vector
                          CVEC Av1(dddm.skeletonsnippetnumber_local_to_global[ssn].size());
                          snippetA.mv(v1,Av1);

                          // now multiply with the rows
                          for (auto& rowentry1 : subdomain_local_to_snippet_local) // loop over all subdomains intersected by the volume snippet
                            for (auto& rowentry2 : rowentry1.second)
                              {
                                auto prow = rowentry1.first;
                                auto krow = rowentry2.first;
                                auto iprojected = rowentry2.second;

                                // restrict basis vector [prow][krow] to snippet
                                CVEC v2(dddm.skeletonsnippetnumber_local_to_global[ssn].size());
                                v2 = 0.0;
                                for (size_t i=0; i<dddm.skeletonsnippetnumber_local_to_global[ssn].size(); i++) // loop over all dofs in snippet
                                  {
                                    auto iglobal = dddm.skeletonsnippetnumber_local_to_global[ssn][i]; // global dof
                                    auto it = dddm.subdomainnumber_global_to_local[iglobal].find(prow);
                                    if (it!=dddm.subdomainnumber_global_to_local[iglobal].end())
                                      v2[i] = coarse_subdomain_basis_vectors[targetlevel][prow][krow][it->second]; // pick entry from basis vector
                                  }

                                // now we have one entry
                                projectedA[iprojected][jprojected] = v2 * Av1;
                              }
                        }
                    if (verbose_>1) std::cout << " ... done" << std::endl;
                  } // loop over skeleton snippets
              } // if there are skeleton snippets
          } // if we need to project snippet matrices
     }

    // allocated defect and correction vectors for the solver
    dcoarse.resize(Number_Of_Levels);
    vcoarse.resize(Number_Of_Levels);
    for (int l=1; l<Number_Of_Levels; l++)
      {
        dcoarse[l] = std::shared_ptr<CVEC>(new CVEC(dddmhierarchy[l]->n_global_dofs));
        vcoarse[l] = std::shared_ptr<CVEC>(new CVEC(dddmhierarchy[l]->n_global_dofs));
      }

    setup_time_total = timer_total.elapsed();

    // print info
    for (size_t l=0; l<Number_Of_Levels; l++)
      {
        size_t sum=0;
        size_t maximum=0;
        for (size_t p=0; p<nsubdomains[l]; p++)
          {
            sum += dddmhierarchy[l]->subdomainnumber_local_to_global[p].size();
            maximum = std::max(maximum,dddmhierarchy[l]->subdomainnumber_local_to_global[p].size());
          }
        std::cout << "LEVEL " << l << ": " << dddmhierarchy[l]->subdomainnumber_global_to_local.size() << " dofs in " << nsubdomains[l] << " subdomains";
        std::cout << " max=" << maximum;
        std::cout << " avg=" << ((double)sum)/nsubdomains[l];
        std::cout << " ideal=" << ((double)dddmhierarchy[l]->subdomainnumber_global_to_local.size())/nsubdomains[l];
        std::cout << std::endl;
      }

    // determine unaccounted time
    double sum_of_parts = setup_time_sequential;
    for (auto& t : setup_time_parallelizable)
      sum_of_parts += t;
    for (auto& v : setup_time_procs)
      for (auto& t : v)
        sum_of_parts += t;
    std::cout << "TIMING SETUP PHASE total=" << setup_time_total << " sum of parts=" << sum_of_parts << " unaccounted=" << setup_time_total-sum_of_parts << " %="
              << 100.0*(setup_time_total-sum_of_parts)/setup_time_total << std::endl;
    std::cout << "TIMING SETUP PHASE sequential part=" << setup_time_sequential << std::endl;

    // compute estimated parallel runtime
    double estimated_runtime = setup_time_sequential;
    for (size_t l=0; l<Number_Of_Levels; l++)
      {
        estimated_runtime += setup_time_parallelizable[l]/nsubdomains[l];
        double maxtime=0.0;
        for (auto& t : setup_time_procs[l])
          maxtime = std::max(maxtime,t);
        double mintime=1e100;
        for (auto& t : setup_time_procs[l])
          mintime = std::min(mintime,t);
        estimated_runtime += maxtime;
        std::cout << "TIMING SETUP PHASE parallel part level=" << l << " parallelizable=" << setup_time_parallelizable[l]/nsubdomains[l]
                  << " min=" << mintime
                  << " max=" << maxtime
                  << std::endl;
      }
    std::cout << "TIMING SETUP PHASE Tseq=" << sum_of_parts << " Tpar=" << estimated_runtime << " S=" << sum_of_parts/estimated_runtime << std::endl;

    double coarse_umfpack_time = 0.0;
    for (size_t l=0; l<Number_Of_Levels; l++)
      {
        double maxtime_umf=0.0;
        for (auto& t : setup_time_umfpack[l])
          maxtime_umf = std::max(maxtime_umf,t);
        if (l>0) coarse_umfpack_time += maxtime_umf;
        double maxtime_par=0.0;
        for (auto& t : setup_time_procs[l])
          maxtime_par = std::max(maxtime_par,t);
        std::cout << "TIMING SETUP PHASE umfpack vs other level=" << l << " ev+umf=" << maxtime_par << " umf=" << maxtime_umf << std::endl;
      }
    std::cout << "TIMING SETUP PHASE total coarse umfpack time=" << coarse_umfpack_time << std::endl;

    // timings for coarse eigensolves
    for (size_t l=0; l<Number_Of_Levels; l++)
      {
        double setup_gevp=0.0;
        double QR=0.0;
        double projection=0.0;
        double gevp=0.0;
        for (auto& t : setup_time_setup_gevp[l]) setup_gevp =  std::max(t,setup_gevp);
        for (auto& t : setup_time_QR[l]) QR =  std::max(t,QR);
        for (auto& t : setup_time_projection[l]) projection =  std::max(t,projection);
        for (auto& t : setup_time_gevp[l]) gevp =  std::max(t,gevp);
        std::cout << "TIMING SETUP PHASE coarse GEVP level=" << l
                  << " setup=" << setup_gevp
                  << " QR=" << QR
                  << " projection=" << projection
                  << " gevp=" << gevp
                  << std::endl;
      }
  } // end of constructor

  //! \brief get number of basis vectors in subdomain
  size_t number_of_basis_vectors (int subdomain)
  {
    if (subdomain<0 || subdomain>=nsubdomains[0])
      DUNE_THROW(Dune::Exception, "number_of_basis_vectors: subdomain out of range.");
    return fine_subdomain_basis_vectors[subdomain].size();
  }

  //! \brief get number of basis vectors in subdomain
  size_t number_of_basis_vectors (int level, int subdomain)
  {
    if (level<0 || level>=Number_Of_Levels-1)
      DUNE_THROW(Dune::Exception, "number_of_basis_vectors: level out of range.");
    if (subdomain<0 || subdomain>=nsubdomains[level])
      DUNE_THROW(Dune::Exception, "number_of_basis_vectors: subdomain out of range.");
    if (level==0)
      return fine_subdomain_basis_vectors[subdomain].size();
    else
      return coarse_subdomain_basis_vectors[level][subdomain].size();
  }

  //! \brief get number of basis vectors in subdomain
  size_t number_of_kerAB_vectors (int level, int subdomain)
  {
    if (level<1 || level>=Number_Of_Levels-1)
      DUNE_THROW(Dune::Exception, "number_of_basis_vectors: level out of range.");
    if (subdomain<0 || subdomain>=nsubdomains[level])
      DUNE_THROW(Dune::Exception, "number_of_basis_vectors: subdomain out of range.");
    return coarse_kerAB_vectors[level][subdomain].size();
  }

  //! get number of levels
  size_t number_of_levels ()
  {
    return Number_Of_Levels;
  }

    //! get number of levels
  size_t number_of_subdomains (int level)
  {
    if (level<0 || level>=Number_Of_Levels-1)
      DUNE_THROW(Dune::Exception, "number_of_basis_vectors: level out of range.");
    return nsubdomains[level];
  }

  //! get basis vector k in subdomain p
  VEC prolongated_basis_vector (int p, int k)
  {
    if (p<0 && p>=nsubdomains[0])
      DUNE_THROW(Dune::Exception, "number_of_basis_vectors: subdomain out of range.");

    VEC v(dddmhierarchy[0]->subdomainnumber_global_to_local.size()); // the global vector to be returned
    v = 0.0;
    auto n = dddmhierarchy[0]->subdomainnumber_local_to_global[p].size();
    for (size_t i=0; i<n; i++)
      v[dddmhierarchy[0]->subdomainnumber_local_to_global[p][i]] = fine_subdomain_basis_vectors[p][k][i];
    return v;
  }

  //! get partition of unity in subdomain p
  VEC prolongated_partition_of_unity (int p)
  {
    if (p<0 && p>=nsubdomains[0])
      DUNE_THROW(Dune::Exception, "number_of_basis_vectors: subdomain out of range.");

    VEC v(dddmhierarchy[0]->subdomainnumber_global_to_local.size()); // the global vector to be returned
    v = 0.0;
    auto n = dddmhierarchy[0]->subdomainnumber_local_to_global[p].size();
    for (size_t i=0; i<n; i++)
      v[dddmhierarchy[0]->subdomainnumber_local_to_global[p][i]] = partition_of_unity[0][p][i];
    return v;
  }

  //! get basis vector k in subdomain p on level l
  VEC prolongated_basis_vector (int l, int p, int k)
  {
    if (l<0 || l>=Number_Of_Levels-1)
      DUNE_THROW(Dune::Exception, "prolongated_basis_vector: level out of range.");
    if (p<0 || p>=nsubdomains[l])
      DUNE_THROW(Dune::Exception, "prolongated_basis_vector subdomain out of range.");
    if (k<0 && k>=coarse_subdomain_basis_vectors[l][p].size())
      DUNE_THROW(Dune::Exception, "prolongated_basis_vector: k out of range.");

    // make the result
    VEC result(dddmhierarchy[0]->subdomainnumber_global_to_local.size()); // the global vector to be returned
    result = 0.0;

    if (l==0)
      {
        // we are already on the finest level. No prolongation necessary
        auto n = dddmhierarchy[0]->subdomainnumber_local_to_global[p].size();
        for (size_t i=0; i<n; i++)
          result[dddmhierarchy[0]->subdomainnumber_local_to_global[p][i]] = fine_subdomain_basis_vectors[p][k][i];
        return result;
      }

    // we are on a coarse level that needs prolongation.

    // make a global dof vector (v) on level l containing the basis vector p,k
    CVEC v(dddmhierarchy[l]->subdomainnumber_global_to_local.size()); // the global vector to be prolonged
    v = 0.0;
    {
      auto nsubdomain = dddmhierarchy[l]->subdomainnumber_local_to_global[p].size();
      for (size_t i=0; i<nsubdomain; i++)
        v[dddmhierarchy[l]->subdomainnumber_local_to_global[p][i]] = coarse_subdomain_basis_vectors[l][p][k][i];
    }

    for (int targetlevel=l-1; targetlevel>=0; targetlevel--)
      if (targetlevel>0)
        {
          // intermediate level
          CVEC w(dddmhierarchy[targetlevel]->subdomainnumber_global_to_local.size()); // intermediate vector on this level
          w = 0.0;

          // prolongate one level
          for (size_t i=0; i<v.size(); i++) // loop over all dofs from targetlevel+1
            {
              auto p = dddmhierarchy[targetlevel]->basisvector_global_to_subdomainnumber[i];
              auto k = dddmhierarchy[targetlevel]->basisvector_global_to_local[i];
              // add basis vector to result
              auto n = dddmhierarchy[targetlevel]->subdomainnumber_local_to_global[p].size();
              for (size_t j=0; j<n; j++)
                w[dddmhierarchy[targetlevel]->subdomainnumber_local_to_global[p][j]] += coarse_subdomain_basis_vectors[targetlevel][p][k][j]*v[i];
            }

          // w is the input on the next level
          v = w;
        }
      else
        {
          // targetlevel==0; prolongate to result
          for (size_t i=0; i<v.size(); i++) // loop over all dofs from targetlevel+1
            {
              auto p = dddmhierarchy[targetlevel]->basisvector_global_to_subdomainnumber[i];
              auto k = dddmhierarchy[targetlevel]->basisvector_global_to_local[i];
              // add basis vector to result
              auto n = dddmhierarchy[targetlevel]->subdomainnumber_local_to_global[p].size();
              for (size_t j=0; j<n; j++)
                {
                  auto blub = fine_subdomain_basis_vectors[p][k][j];
                  blub *= v[i][0];
                  result[dddmhierarchy[targetlevel]->subdomainnumber_local_to_global[p][j]] += blub;
                }
            }
        }

    return result;
  }

    //! get basis vector k in subdomain p on level l
  VEC prolongated_kerAB_vector (int l, int p, int k)
  {
    if (l<1 || l>=Number_Of_Levels-1)
      DUNE_THROW(Dune::Exception, "prolongated_kerAB_vector: level out of range.");
    if (p<0 || p>=nsubdomains[l])
      DUNE_THROW(Dune::Exception, "prolongated_kerAB_vector subdomain out of range.");
    if (k<0 && k>=coarse_kerAB_vectors[l][p].size())
      DUNE_THROW(Dune::Exception, "prolongated_kerAB_vector: k out of range.");

    // make the result
    VEC result(dddmhierarchy[0]->subdomainnumber_global_to_local.size()); // the global vector to be returned
    result = 0.0;

    // make a global dof vector (v) on level l containing the basis vector p,k
    CVEC v(dddmhierarchy[l]->subdomainnumber_global_to_local.size()); // the global vector to be prolonged
    v = 0.0;
    {
      auto nsubdomain = dddmhierarchy[l]->subdomainnumber_local_to_global[p].size();
      for (size_t i=0; i<nsubdomain; i++)
        v[dddmhierarchy[l]->subdomainnumber_local_to_global[p][i]] = coarse_kerAB_vectors[l][p][k][i];
    }

    for (int targetlevel=l-1; targetlevel>=0; targetlevel--)
      if (targetlevel>0)
        {
          // intermediate level
          CVEC w(dddmhierarchy[targetlevel]->subdomainnumber_global_to_local.size()); // intermediate vector on this level
          w = 0.0;

          // prolongate one level
          for (size_t i=0; i<v.size(); i++) // loop over all dofs from targetlevel+1
            {
              auto p = dddmhierarchy[targetlevel]->basisvector_global_to_subdomainnumber[i];
              auto k = dddmhierarchy[targetlevel]->basisvector_global_to_local[i];
              // add basis vector to result
              auto n = dddmhierarchy[targetlevel]->subdomainnumber_local_to_global[p].size();
              for (size_t j=0; j<n; j++)
                w[dddmhierarchy[targetlevel]->subdomainnumber_local_to_global[p][j]] += coarse_subdomain_basis_vectors[targetlevel][p][k][j]*v[i];
            }

          // w is the input on the next level
          v = w;
        }
      else
        {
          // targetlevel==0; prolongate to result
          for (size_t i=0; i<v.size(); i++) // loop over all dofs from targetlevel+1
            {
              auto p = dddmhierarchy[targetlevel]->basisvector_global_to_subdomainnumber[i];
              auto k = dddmhierarchy[targetlevel]->basisvector_global_to_local[i];
              // add basis vector to result
              auto n = dddmhierarchy[targetlevel]->subdomainnumber_local_to_global[p].size();
              for (size_t j=0; j<n; j++)
                {
                  auto blub = fine_subdomain_basis_vectors[p][k][j];
                  blub *= v[i][0];
                  result[dddmhierarchy[targetlevel]->subdomainnumber_local_to_global[p][j]] += blub;
                }
            }
        }

    return result;
  }

  //! \brief export information from domain decomposition
  const std::vector<std::set<size_t>>& get_overlappingsubdomains (int l)
  {
    return ddhierarchy[l]->overlappingsubdomains;
  }

  //! \brief export information from domain decomposition
  const std::vector<size_t>& get_element_to_volumesnippetnumber (int l)
  {
    return ddhierarchy[l]->element_to_volumesnippetnumber;
  }

  //! \brief export information from domain decomposition
  const std::vector<size_t>& get_partition (int l)
  {
    return ddhierarchy[l]->partition;
  }

  //! \brief export information from domain decomposition
  const std::vector<size_t>& get_vizpartition (int l)
  {
    return ddhierarchy[l]->vizpartition;
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
    if (cycle=="mras") // multiplicative over levels, RAS within levels
      {
        // we start with the fine level subdomain solves
        for (size_t p=0; p<nsubdomains[0]; p++)
          {
            // dof blocks in subdomain p
            auto n = dddmhierarchy[0]->subdomainnumber_local_to_global[p].size();

            // set up right hand side in subdomain
            VEC dlocal(n);
            for (size_t i=0; i<n; i++)
              {
                dlocal[i] = d[dddmhierarchy[0]->subdomainnumber_local_to_global[p][i]];
                dlocal[i] *= subdomain_boundary_mask[0][p][i];
              }

            // solve subdomain problem
            Dune::InverseOperatorResult stat;
            VEC vlocal(n);
            vlocal = 0.0;
            finesubdomainsolvers[p]->apply(vlocal,dlocal,stat);

            // multiply with partition of unity -> restricted additive schwarz
            for (size_t i=0; i<n; i++)
              vlocal[i] *= partition_of_unity[0][p][i];

            // accumulate correction from subdomain
            for (size_t i=0; i<n; i++)
              v[dddmhierarchy[0]->subdomainnumber_local_to_global[p][i]] += vlocal[i];
          }
        // now the fine subdomains are completed

        // update defect
        VEC d2(d);
        pAglobal->mmv(v,d2);

        // restrict defect to level 1; defect vectors have been allocated already
        // due to different types, restriction level 0->1 has to be handled seperately
        {
          *dcoarse[1] = 0.0; // clear defect on level 1
          DomainDecompositionDOFMapper& finedddm = *dddmhierarchy[0];
          for (size_t p=0; p<nsubdomains[0]; p++) // loop over subdomains one level finer
            for (size_t k=0; k<finedddm.basisvector_local_to_global[p].size(); k++) // loop over all basis vectors
              {
                auto icoarse = finedddm.basisvector_local_to_global[p][k]; // associated global degree of freedom on the coarser level
                auto n = finedddm.subdomainnumber_local_to_global[p].size(); // size of the basis vector
                for (size_t i=0; i<n; i++)
                  (*dcoarse[1])[icoarse] += d2[finedddm.subdomainnumber_local_to_global[p][i]]*fine_subdomain_basis_vectors[p][k][i];
              }
        }

        // perform downward part of V-cycle
        for (int targetlevel=1; targetlevel<Number_Of_Levels-1; targetlevel++)
          {
            // info for that level
            DomainDecompositionDOFMapper& dddm = *dddmhierarchy[targetlevel];

            // clear global correction on this level
            *vcoarse[targetlevel] = 0.0;

            // and perform subdomain solves
            for (size_t p=0; p<nsubdomains[targetlevel]; p++)
              {
                // solve in subdomain p
                auto n = dddm.subdomainnumber_local_to_global[p].size();

                // set up right hand side in subdomain (pick out)
                CVEC dlocal(n);
                for (size_t i=0; i<n; i++)
                  {
                    dlocal[i] = (*dcoarse[targetlevel])[dddm.subdomainnumber_local_to_global[p][i]]; // this is a picking out operator
                    dlocal[i] *= subdomain_boundary_mask[targetlevel][p][i]; // the Dirichlet condition
                  }

                // solve subdomain problem
                Dune::InverseOperatorResult stat;
                CVEC vlocal(n);
                vlocal = 0.0;
                coarsesubdomainsolvers[targetlevel][p]->apply(vlocal,dlocal,stat);

                // multiply with partition of unity -> restricted additive schwarz
                for (size_t i=0; i<n; i++)
                  vlocal[i] *= partition_of_unity[targetlevel][p][i];

                // accumulate correction from subdomain to global correction on this level
                for (size_t i=0; i<n; i++)
                  (*vcoarse[targetlevel])[dddm.subdomainnumber_local_to_global[p][i]] += vlocal[i]; // this is the extension operator
              }

              // update defect on that level
              global_coarse_matrix[targetlevel]->mmv(*vcoarse[targetlevel],*dcoarse[targetlevel]);

              // restrict defect to next coarser level (which must exist)
              *dcoarse[targetlevel+1] = 0.0; // clear defect on targetlevel
              for (size_t p=0; p<nsubdomains[targetlevel]; p++) // loop over subdomains
                for (size_t k=0; k<dddm.basisvector_local_to_global[p].size(); k++) // loop over all basis vectors in that subdomain
                  {
                    auto icoarse = dddm.basisvector_local_to_global[p][k]; // associated global degree of freedom on the coarser level
                    auto n = dddm.subdomainnumber_local_to_global[p].size(); // size of the basis vector
                    for (size_t i=0; i<n; i++)
                      (*dcoarse[targetlevel+1])[icoarse] += (*dcoarse[targetlevel])[dddm.subdomainnumber_local_to_global[p][i]]*coarse_subdomain_basis_vectors[targetlevel][p][k][i];
                  }
          }

        // coarsest level solve
        {
          // the coarsest level; there is only one subdomain
          auto targetlevel = Number_Of_Levels-1;
          Dune::InverseOperatorResult stat;
          *vcoarse[targetlevel] = 0.0;
          globalcoarsesolvers[targetlevel]->apply(*vcoarse[targetlevel],*dcoarse[targetlevel],stat);
        }

        // upward part of V-cycle
        for (int sourcelevel=Number_Of_Levels-1; sourcelevel>1; sourcelevel--) // loop over remaining levels
          {
            DomainDecompositionDOFMapper& finedddm = *dddmhierarchy[sourcelevel-1]; // targetlevel
            for (size_t p=0; p<nsubdomains[sourcelevel-1]; p++) // loop over subdomains one level finer
              for (size_t k=0; k<finedddm.basisvector_local_to_global[p].size(); k++) // loop over all basis vectors
                {
                  auto icoarse = finedddm.basisvector_local_to_global[p][k]; // associated global degree of freedom on the coarser level
                  auto n = finedddm.subdomainnumber_local_to_global[p].size(); // size of the basis vector
                  for (size_t i=0; i<n; i++)
                    (*vcoarse[sourcelevel-1])[finedddm.subdomainnumber_local_to_global[p][i]] += (*vcoarse[sourcelevel])[icoarse][0]*coarse_subdomain_basis_vectors[sourcelevel-1][p][k][i];
                }
          }
        // prolongation 1(=sourcelevel)->0
        {
          DomainDecompositionDOFMapper& finedddm = *dddmhierarchy[0]; // targetlevel
          for (size_t p=0; p<nsubdomains[0]; p++) // loop over subdomains one level finer
            for (size_t k=0; k<finedddm.basisvector_local_to_global[p].size(); k++) // loop over all basis vectors
              {
                auto icoarse = finedddm.basisvector_local_to_global[p][k]; // associated global degree of freedom on the coarser level
                auto n = finedddm.subdomainnumber_local_to_global[p].size(); // size of the basis vector
                for (size_t i=0; i<n; i++)
                  v[finedddm.subdomainnumber_local_to_global[p][i]] += (*vcoarse[1])[icoarse][0]*fine_subdomain_basis_vectors[p][k][i];
              }
        }

        return; // end of mras cycle
      }

    //===========================================
    // now the default cycle form: fully additive
    //===========================================

    // we start with the fine level subdomain solves
    for (size_t p=0; p<nsubdomains[0]; p++)
      {
        // dof blocks in subdomain p
        auto n = dddmhierarchy[0]->subdomainnumber_local_to_global[p].size();

        // set up right hand side in subdomain
        VEC dlocal(n);
        for (size_t i=0; i<n; i++)
          {
            dlocal[i] = d[dddmhierarchy[0]->subdomainnumber_local_to_global[p][i]];
            dlocal[i] *= subdomain_boundary_mask[0][p][i];
          }

        // solve subdomain problem
        Dune::InverseOperatorResult stat;
        VEC vlocal(n);
        vlocal = 0.0;
        finesubdomainsolvers[p]->apply(vlocal,dlocal,stat);

        // multiply with partition of unity -> restricted additive schwarz
        if (cycle=="aras")
          for (size_t i=0; i<n; i++)
            vlocal[i] *= partition_of_unity[0][p][i];

        // accumulate correction from subdomain
        for (size_t i=0; i<n; i++)
          v[dddmhierarchy[0]->subdomainnumber_local_to_global[p][i]] += vlocal[i];
      }
    // now the fine subdomains are completed

    // restrict defect to all coarse levels; defect vectors have been allocated already
    // due to different types, restriction level 0->1 has to be handled seperately
    {
      *dcoarse[1] = 0.0; // clear defect on level 1
      DomainDecompositionDOFMapper& finedddm = *dddmhierarchy[0];
      for (size_t p=0; p<nsubdomains[0]; p++) // loop over subdomains one level finer
        for (size_t k=0; k<finedddm.basisvector_local_to_global[p].size(); k++) // loop over all basis vectors
          {
            auto icoarse = finedddm.basisvector_local_to_global[p][k]; // associated global degree of freedom on the coarser level
            auto n = finedddm.subdomainnumber_local_to_global[p].size(); // size of the basis vector
            for (size_t i=0; i<n; i++)
              (*dcoarse[1])[icoarse] += d[finedddm.subdomainnumber_local_to_global[p][i]]*fine_subdomain_basis_vectors[p][k][i];
          }
    }
    for (int targetlevel=2; targetlevel<Number_Of_Levels; targetlevel++) // loop over remaining levels
      {
        *dcoarse[targetlevel] = 0.0; // clear defect on targetlevel
        DomainDecompositionDOFMapper& finedddm = *dddmhierarchy[targetlevel-1];
        for (size_t p=0; p<nsubdomains[targetlevel-1]; p++) // loop over subdomains one level finer
          for (size_t k=0; k<finedddm.basisvector_local_to_global[p].size(); k++) // loop over all basis vectors
            {
              auto icoarse = finedddm.basisvector_local_to_global[p][k]; // associated global degree of freedom on the coarser level
              auto n = finedddm.subdomainnumber_local_to_global[p].size(); // size of the basis vector
              for (size_t i=0; i<n; i++)
                (*dcoarse[targetlevel])[icoarse] += (*dcoarse[targetlevel-1])[finedddm.subdomainnumber_local_to_global[p][i]]*coarse_subdomain_basis_vectors[targetlevel-1][p][k][i];
            }
      }

    // solve subdomain problems
    for (int targetlevel=1; targetlevel<Number_Of_Levels; targetlevel++) // loop over remaining levels
      if (targetlevel==Number_Of_Levels-1)
        {
          // the coarsest level
          Dune::InverseOperatorResult stat;
          *vcoarse[targetlevel] = 0.0;
          globalcoarsesolvers[targetlevel]->apply(*vcoarse[targetlevel],*dcoarse[targetlevel],stat);
        }
      else
        {
          // clear global correection on this level
          *vcoarse[targetlevel] = 0.0;

          // an intermediate level
          for (size_t p=0; p<nsubdomains[targetlevel]; p++)
            {
              // solve in subdomain p
              DomainDecompositionDOFMapper& dddm = *dddmhierarchy[targetlevel];
              auto n = dddm.subdomainnumber_local_to_global[p].size();

              // set up right hand side in subdomain (pick out)
              CVEC dlocal(n);
              for (size_t i=0; i<n; i++)
                {
                  dlocal[i] = (*dcoarse[targetlevel])[dddm.subdomainnumber_local_to_global[p][i]]; // this is a picking out operator
                  dlocal[i] *= subdomain_boundary_mask[targetlevel][p][i]; // the Dirichlet condition
                }

              // solve subdomain problem
              Dune::InverseOperatorResult stat;
              CVEC vlocal(n);
              vlocal = 0.0;
              coarsesubdomainsolvers[targetlevel][p]->apply(vlocal,dlocal,stat);

              // multiply with partition of unity -> restricted additive schwarz
              if (cycle=="aras")
                for (size_t i=0; i<n; i++)
                  vlocal[i] *= partition_of_unity[targetlevel][p][i];

              // accumulate correction from subdomain to global correction on this level
              for (size_t i=0; i<n; i++)
                (*vcoarse[targetlevel])[dddm.subdomainnumber_local_to_global[p][i]] += vlocal[i]; // this is the extension operator
            }
        }

    // now we have all the correction on the coarse levels in global vectors
    // prolongation
    for (int sourcelevel=Number_Of_Levels-1; sourcelevel>1; sourcelevel--) // loop over remaining levels
      {
        DomainDecompositionDOFMapper& finedddm = *dddmhierarchy[sourcelevel-1]; // targetlevel
        for (size_t p=0; p<nsubdomains[sourcelevel-1]; p++) // loop over subdomains one level finer
          for (size_t k=0; k<finedddm.basisvector_local_to_global[p].size(); k++) // loop over all basis vectors
            {
              auto icoarse = finedddm.basisvector_local_to_global[p][k]; // associated global degree of freedom on the coarser level
              auto n = finedddm.subdomainnumber_local_to_global[p].size(); // size of the basis vector
              for (size_t i=0; i<n; i++)
                (*vcoarse[sourcelevel-1])[finedddm.subdomainnumber_local_to_global[p][i]] += (*vcoarse[sourcelevel])[icoarse][0]*coarse_subdomain_basis_vectors[sourcelevel-1][p][k][i];
            }
      }
    // prolongation 1(=sourcelevel)->0
    {
        DomainDecompositionDOFMapper& finedddm = *dddmhierarchy[0]; // targetlevel
        for (size_t p=0; p<nsubdomains[0]; p++) // loop over subdomains one level finer
          for (size_t k=0; k<finedddm.basisvector_local_to_global[p].size(); k++) // loop over all basis vectors
            {
              auto icoarse = finedddm.basisvector_local_to_global[p][k]; // associated global degree of freedom on the coarser level
              auto n = finedddm.subdomainnumber_local_to_global[p].size(); // size of the basis vector
              for (size_t i=0; i<n; i++)
                v[finedddm.subdomainnumber_local_to_global[p][i]] += (*vcoarse[1])[icoarse][0]*fine_subdomain_basis_vectors[p][k][i];
            }
    }
  } // end of apply

  /*!
    \brief Clean up.
  */
  virtual void post (VEC& x)
  {
  }

};


#endif // HAVE_SUITESPARSE_UMFPACK || DOXYGEN

#endif // Udune_ftworth_HH
