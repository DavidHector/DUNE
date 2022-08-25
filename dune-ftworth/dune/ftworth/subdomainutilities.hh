#ifndef Udune_ftworth_subdomainutilities_HH
#define Udune_ftworth_subdomainutilities_HH

#include<vector>
#include<map>
#include<unordered_map>
#include<set>

/** \brief Provides information about the degrees of freedom in the subdomain on the finest level
 *
 * \author Peter Bastian
 *
 * \tparam GFS grid function space
 *
 * Collects all information needed to assemble the subdomain problems
 */
template<typename GFS>
class FineLevelSubdomainInformation
{
public:
  // types
  using GV = typename GFS::Traits::GridView;
  using IS = typename GV::IndexSet;
  using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
  using LFSCache = Dune::PDELab::LFSIndexCache<LFS>;
  using RF = size_t;
  using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
  using ISTLV = Dune::PDELab::Backend::Native<Z>;
  using VIEW = typename Z::template LocalView<LFSCache>;

  // the information provided by this class
  const GV& gv; // from the gfs
  const IS& indexset; // from the gv
  LFS lfs; // a local function space for the gfs
  LFSCache lfscache; // lfs cache for accessing dofs of an element
  Z countindices; // one-dimensional flat index for each dof
  VIEW view; // local part of the vector associated with an element
  unsigned int subdomains; // the number of subdomains
  std::shared_ptr<std::vector<std::vector<size_t>>> local_to_global; // local_to_global[p][i] gives global index for local index i in subdomain p
  std::shared_ptr<std::vector<std::set<size_t>>> global_to_subdomains; // global_to_subdomains[i] gives a set of subdomain numbers for global index i
  std::shared_ptr<std::vector<std::map<size_t,size_t>>> global_to_local; // global_to_local[i][p] gives the local index in subdomain p for global index i (make sure that subdomain p contains i!)
  RF blocksize; // assume that vector consists of blocks of this fixed size

  FineLevelSubdomainInformation (const GFS& gfs, unsigned partitions, std::vector<std::set<size_t>> overlappingsubdomains)
    : gv(gfs.gridView()), indexset(gv.indexSet()), lfs(gfs), lfscache(lfs), countindices(gfs), view(countindices), subdomains(partitions)
  {
    // fill z with consecutive indices
    // we assume that our ISTL vector is a BlockVector<FieldVector<>>
    ISTLV& istlz = Dune::PDELab::Backend::native(countindices);
    RF n = 0; 
    for (int i=0; i<istlz.size(); i++)
      for (int j=0; j<istlz[i].size(); j++)
        istlz[i][j] = n++;

    // determine block size
    blocksize = istlz[0].size();

    // indices in each subdomain
    std::vector<std::set<size_t>> indices_per_subdomain(partitions); // maps partition to a set of global indices (not exported)
    for (const auto& e : elements(gv))
      {
        lfs.bind(e);
        lfscache.update();
        std::vector<RF> localdata(lfs.size());
        view.bind(lfscache);
        view.read(localdata);
        view.unbind();
        for (const auto& p : overlappingsubdomains[indexset.index(e)])
          for (size_t i = 0; i<lfs.size(); ++i)
            indices_per_subdomain[p].insert(localdata[i]/blocksize);
      }

    // map local to global indices
    local_to_global = std::shared_ptr<std::vector<std::vector<size_t>>>(new std::vector<std::vector<size_t>>(partitions));
    for (int p=0; p<partitions; p++)
      for (auto& i : indices_per_subdomain[p])
        (*local_to_global)[p].push_back(i);
    // for (int p=0; p<partitions; p++)
    //   std::cout << "subdomain " << p << " has " << (*local_to_global)[p].size() << " dof blocks" << std::endl;

    // for each global index store the partitions containing this index
    // this allows us also to set up the subdomain matrices
    global_to_subdomains = std::shared_ptr<std::vector<std::set<size_t>>>(new std::vector<std::set<size_t>>(istlz.size()));
    for (const auto& e : elements(gv))
      {
        lfs.bind(e);
        lfscache.update();
        std::vector<RF> localdata(lfs.size());
        view.bind(lfscache);
        view.read(localdata);
        view.unbind();
        for (size_t i = 0; i<lfs.size(); ++i)
          {
            int globalindex = localdata[i]/blocksize;
            for (const auto& p : overlappingsubdomains[indexset.index(e)])
              (*global_to_subdomains)[globalindex].insert(p);
          }
      }

    // map global to local indices
    // you have to check before whether the global index is used in that subdomain
    global_to_local = std::shared_ptr<std::vector<std::map<size_t,size_t>>>(new std::vector<std::map<size_t,size_t>>(istlz.size()));
    for (size_t p=0; p<partitions; ++p)
      for (size_t ilocal=0; ilocal<(*local_to_global)[p].size(); ++ilocal)
        (*global_to_local)[(*local_to_global)[p][ilocal]][p] = ilocal;
  }
};

// set up subdomain matrices
template<typename ISTLM>
std::vector<std::shared_ptr<ISTLM>> set_up_subdomain_matrices (const ISTLM& istlA,
							       const std::vector<std::vector<size_t>>& local_to_global,
							       const std::vector<std::set<size_t>>& global_to_subdomains,
							       const std::vector<std::map<size_t,size_t>>& global_to_local)
{
  // build up sparse matrices
  auto subdomains = local_to_global.size();
  int avg = istlA.nonzeroes()/istlA.N()+2;
  std::vector<std::shared_ptr<ISTLM>> matrices;
  for (size_t p=0; p<subdomains; ++p)
    {
      auto n = local_to_global[p].size(); // number of DOFs in subdomain p
      //std::cout << "building matrix in subdomain " << p << " size " << n << " times " << n << " avg row size " << avg << std::endl;
      auto pA = std::shared_ptr<ISTLM>(new ISTLM(n,n,avg,0.15,ISTLM::implicit)); // a new sparse matrix
      for (size_t i=0; i<n; i++)
        {
          auto iglobal = local_to_global[p][i];
          auto cIt = istlA[iglobal].begin();
          auto cEndIt = istlA[iglobal].end();
          for (; cIt!=cEndIt; ++cIt)
            {
              auto jglobal = cIt.index();
              if (global_to_subdomains[jglobal].count(p)==0) continue;
              // ok, now we have an entry
              auto it = global_to_local[jglobal].find(p); // we must find it; this is just because operator[] is non const
              //auto j = global_to_local[jglobal][p]; // global to local
              auto j = it->second; // global to local
              pA->entry(i,j) = 0.0;
            }
        }
      auto stats = pA->compress();
      matrices.push_back(pA);
    }
  return matrices;
}

// implement Dirichlet boundary conditions in subdomain matrices
// returns floating subdomain flag
template<typename ISTLV, typename ISTLM>
std::vector<bool> assemble_dirichlet_constraints (const std::vector<std::vector<size_t>>& local_to_global,
						  const ISTLV& dirichlet,
						  std::vector<std::shared_ptr<ISTLM>>& matrices)
{
  // block size and subdomains
  using BlockMatrix = typename ISTLM::block_type;
  const size_t blocksize = BlockMatrix::rows;
  auto subdomains = local_to_global.size();

  // loop over all matrices
  std::vector<bool> floating(subdomains);
  for (size_t p=0; p<subdomains; ++p)
    {
      ISTLM& A = *matrices[p]; // matrix in subdomain p
      floating[p] = true;
      for (size_t ilocal=0; ilocal<A.N(); ilocal++)
        {
          auto iglobal = local_to_global[p][ilocal];
          for (size_t compi=0; compi<blocksize; compi++)
            {
              //std::cout << "Index " << iglobal << "," << compi << " dirchlet=" << dirichlet[iglobal][compi] << std::endl;
              if (dirichlet[iglobal][compi]==0.0)
                {
                  // we have a dirichlet row
                  floating[p] = false;
                  auto cIt = A[ilocal].begin();
                  auto cEndIt = A[ilocal].end();
                  for (; cIt!=cEndIt; ++cIt)
                    if (cIt.index()==ilocal)
                      for (size_t compj=0; compj<blocksize; compj++) (*cIt)[compi][compj] = (compj==compi) ? 1.0 : 0.0;
                    else
                      for (size_t compj=0; compj<blocksize; compj++) (*cIt)[compi][compj] = 0.0;
                }
              else // this is now for symmetric elimination; the matrix is then block diagonal
                {
                  // regular row
                  auto cIt = A[ilocal].begin();
                  auto cEndIt = A[ilocal].end();
                  for (; cIt!=cEndIt; ++cIt)
                    {
                      auto jglobal = local_to_global[p][cIt.index()];
                      for (size_t compj=0; compj<blocksize; compj++)
                        if (dirichlet[jglobal][compj]==0.0)
                          (*cIt)[compi][compj] = 0.0;
                    }
                }
            }
        }
    }
  return floating;
}

// assemble global dirichlet constraints in a single subdomain
template<typename ISTLV, typename ISTLM>
bool assemble_dirichlet_constraints_in_subdomain (std::vector<size_t>& local_to_global,
						  const ISTLV& dirichlet,
						  std::shared_ptr<ISTLM> matrix)
{
  // block size and subdomains
  using BlockMatrix = typename ISTLM::block_type;
  const size_t blocksize = BlockMatrix::rows;

  ISTLM& A = *matrix; // subdomain
  bool floating = true;
  for (size_t ilocal=0; ilocal<A.N(); ilocal++)
    {
      auto iglobal = local_to_global[ilocal];
      for (size_t compi=0; compi<blocksize; compi++)
	{
	  //std::cout << "Index " << iglobal << "," << compi << " dirchlet=" << dirichlet[iglobal][compi] << std::endl;
	  if (dirichlet[iglobal][compi]==0.0)
	    {
	      // we have a dirichlet row
	      floating = false;
	      auto cIt = A[ilocal].begin();
	      auto cEndIt = A[ilocal].end();
	      for (; cIt!=cEndIt; ++cIt)
		if (cIt.index()==ilocal)
		  for (size_t compj=0; compj<blocksize; compj++) (*cIt)[compi][compj] = (compj==compi) ? 1.0 : 0.0;
		else
		  for (size_t compj=0; compj<blocksize; compj++) (*cIt)[compi][compj] = 0.0;
	    }
	  else // this is now for symmetric elimination; the matrix is then block diagonal
	    {
	      // regular row
	      auto cIt = A[ilocal].begin();
	      auto cEndIt = A[ilocal].end();
	      for (; cIt!=cEndIt; ++cIt)
		{
		  auto jglobal = local_to_global[cIt.index()];
		  for (size_t compj=0; compj<blocksize; compj++)
		    if (dirichlet[jglobal][compj]==0.0)
		      (*cIt)[compi][compj] = 0.0;
		}
	    }
	}
    }

  return floating;
}

// assemble global dirichlet constraints in a single subdomain
template<typename ISTLV, typename ISTLM>
bool assemble_dirichlet_constraints_in_subdomain (std::vector<size_t>& local_to_global,
						  const ISTLV& dirichlet,
						  ISTLM& A)
{
  // block size and subdomains
  using BlockMatrix = typename ISTLM::block_type;
  const size_t blocksize = BlockMatrix::rows;

  bool floating = true;
  for (size_t ilocal=0; ilocal<A.N(); ilocal++)
    {
      auto iglobal = local_to_global[ilocal];
      for (size_t compi=0; compi<blocksize; compi++)
	{
	  //std::cout << "Index " << iglobal << "," << compi << " dirchlet=" << dirichlet[iglobal][compi] << std::endl;
	  if (dirichlet[iglobal][compi]==0.0)
	    {
	      // we have a dirichlet row
	      floating = false;
	      auto cIt = A[ilocal].begin();
	      auto cEndIt = A[ilocal].end();
	      for (; cIt!=cEndIt; ++cIt)
		if (cIt.index()==ilocal)
		  for (size_t compj=0; compj<blocksize; compj++) (*cIt)[compi][compj] = (compj==compi) ? 1.0 : 0.0;
		else
		  for (size_t compj=0; compj<blocksize; compj++) (*cIt)[compi][compj] = 0.0;
	    }
	  else // this is now for symmetric elimination; the matrix is then block diagonal
	    {
	      // regular row
	      auto cIt = A[ilocal].begin();
	      auto cEndIt = A[ilocal].end();
	      for (; cIt!=cEndIt; ++cIt)
		{
		  auto jglobal = local_to_global[cIt.index()];
		  for (size_t compj=0; compj<blocksize; compj++)
		    if (dirichlet[jglobal][compj]==0.0)
		      (*cIt)[compi][compj] = 0.0;
		}
	    }
	}
    }

  return floating;
}



// find out about dofs (blocks) on the boundary of a subdomain
template<typename ISTLM>
std::vector<std::vector<bool>> get_subdomain_boundary_dofs (const ISTLM& istlA,
							    const std::vector<std::vector<size_t>>& local_to_global,
							    const std::vector<std::set<size_t>>& global_to_subdomains)
{
  // initialize result with false
  auto subdomains = local_to_global.size();
  std::vector<std::vector<bool>> boundary_dofs(subdomains);
  for (size_t p=0; p<subdomains; ++p)
    boundary_dofs[p].resize(local_to_global[p].size(),false);

  // a dof is on the boundary if we do not have the complete row of the matrix
  for (size_t p=0; p<subdomains; ++p)
    for (size_t i=0; i<local_to_global[p].size(); i++)
      {
        auto iglobal = local_to_global[p][i];
        auto cIt = istlA[iglobal].begin();
        auto cEndIt = istlA[iglobal].end();
        for (; cIt!=cEndIt; ++cIt) // run through matrix row
          if (global_to_subdomains[cIt.index()].count(p)==0)
            {
              // there is matrix entry in this row that we do not have in subdomain p 
              boundary_dofs[p][i] = true;
              continue;
            }
      }

  return boundary_dofs;
}

// find out about dofs (blocks) on the boundary of a subdomain
template<typename ISTLM>
std::vector<std::vector<bool>> get_subdomain_boundary_dofs (const ISTLM& istlA, // the global matrix on this level
							    const std::vector<std::vector<size_t>>& local_to_global,
							    const std::vector<std::map<size_t,size_t>>& global_to_local)
{
  // initialize result with false
  auto subdomains = local_to_global.size();
  std::vector<std::vector<bool>> boundary_dofs(subdomains);
  for (size_t p=0; p<subdomains; ++p)
    boundary_dofs[p].resize(local_to_global[p].size(),false);

  // a dof is on the boundary if we do not have the complete row of the matrix
  for (size_t p=0; p<subdomains; ++p)
    for (size_t i=0; i<local_to_global[p].size(); i++)
      {
        auto iglobal = local_to_global[p][i];
        auto cIt = istlA[iglobal].begin();
        auto cEndIt = istlA[iglobal].end();
        for (; cIt!=cEndIt; ++cIt) // run through matrix row
          if (global_to_local[cIt.index()].find(p)==global_to_local[cIt.index()].end())
            {
              // there is matrix entry in this row that we do not have in subdomain p 
              boundary_dofs[p][i] = true;
              continue;
            }
      }

  return boundary_dofs;
}


// a simple class representing a weighted graph
// allows construction and coarsening
class WeightedGraph
{
public:
  std::vector<std::set<size_t>> graph; // subdomain_graph[p] gives the vertex connectivity (including p itself, usually)
  std::vector<size_t> vertex_weights; // one number per vertex
  std::vector<std::map<size_t,size_t>> edge_weights; // one number per offdiagonal entry, diagonal entries may be present but are not used

  // constructor without arguments
  WeightedGraph ()
  {}
    
  // build a weighted graph from global to local information (providing joint degrees of freedoms in subdomains)
  WeightedGraph (size_t subdomains, std::vector<std::map<size_t,size_t>> global_to_local)
  {
    graph.resize(subdomains);
    vertex_weights.resize(subdomains,0);
    edge_weights.resize(subdomains);
    for (size_t j=0; j<global_to_local.size(); j++) // loop over global dof numbers
      for (auto itp = global_to_local[j].begin(); itp!=global_to_local[j].end(); ++itp) // loop over subdomain numbers sharing this degree of freedom
	{
	  vertex_weights[itp->first]++; // increment vertex weight
	  for (auto itq = global_to_local[j].begin(); itq!=global_to_local[j].end(); ++itq) // loop over subdomain numbers sharing this degree of freedom
	    {
	      graph[itp->first].insert(itq->first); // couple the two subdomains because they share a degree of freedom
	      if (itp->first!=itq->first) // edge weights are only calculated for offdiagonal edges
		if (edge_weights[itp->first].find(itq->first)==edge_weights[itp->first].end())
		  edge_weights[itp->first][itq->first] = 1; // new entry
		else
		  edge_weights[itp->first][itq->first] += 1; // increment existing entry
	    }
	}
  }
    
  // produce information to construct coarsening of a graph
  // 1) clustering vertices of the given graph
  // 2) and adding overlap
  std::vector<std::set<size_t>> coarsen_graph (MPI_Comm comm, int parts, int overlap, const std::vector<bool>& disconnected, bool merge_disconnected)
  {
    std::vector<std::set<size_t>> coarse_to_fine_map(parts); // the output
      
    if (parts>1)
      {
	// then we have to call parmetis
	// compute overlapping clusters of vertices of the current graph
	auto graph_partitioning = parmetis_sequential_graph_partitioning(graph,vertex_weights,edge_weights,comm,parts);
	auto overlapping_graph_partitioning = sequential_partition_growing(graph,graph_partitioning,overlap,disconnected,merge_disconnected);

	// at this point overlapping_graph_partitioning[i] gives a set of coarse vertices for each fine vertex
	// now compute the mapping of coarse vertices to a set of fine vertices
	for (size_t finep=0; finep<overlapping_graph_partitioning.size(); finep++) // run over all vertices in fine graph
	  for (auto coarsep : overlapping_graph_partitioning[finep]) // run over all coarse vertices the fine vertex is contained in
	    coarse_to_fine_map[coarsep].insert(finep);
      }
    else
      for (size_t p=0; p<graph.size(); p++)
	coarse_to_fine_map[0].insert(p); // all are in one vertex

    // return result
    return coarse_to_fine_map;
  }

  // given a coarse to fine map describing an overlapping clustering of the given
  // graph, compute a reduced weighted graph
  WeightedGraph reduce (const std::vector<std::set<size_t>>& coarse_to_fine_map)
  {
    // sizes
    auto nf = graph.size();
    auto nc = coarse_to_fine_map.size(); // number of vertices in reduced graph

    // rebuild fine to coarse map (we had this before ...)
    std::vector<std::set<size_t>> fine_to_coarse_map(nf);
    for (size_t p=0; p<nc; p++) // loop over coarse vertices
      for (auto q : coarse_to_fine_map[p]) // loop over all fine vertices contributing to this coarse vertex
	fine_to_coarse_map[q].insert(p);
      
    // now create the reduced graph
    WeightedGraph cg; // the new graph is initially empty
    cg.graph.resize(nc);
    cg.vertex_weights.resize(nc);
    cg.edge_weights.resize(nc);

    for (size_t pfine=0; pfine<nf; ++pfine) // loop over all fine vertices
      {
	// contribution to vertex weight
	for (auto pcoarse : fine_to_coarse_map[pfine]) // loop over all coarse vertices where pfine contributes
	  cg.vertex_weights[pcoarse] += vertex_weights[pfine]; // add fine weight to coarse weight

	// connections
	for (auto qfine : graph[pfine]) // loop over all neighbors of this vertex (on the fine grid)
	  // so pfine and qfine have joint degrees of freedom because they are neighbors
	  for (auto pcoarse : fine_to_coarse_map[pfine]) // pcoarse contains pfine
	    for (auto qcoarse : fine_to_coarse_map[qfine]) // qcoarse contains qfine
	      {
		cg.graph[pcoarse].insert(qcoarse); // so these are neighbors
		// add edge weight only if coarse vertices are different
		if (pcoarse!=qcoarse)
		  {
		    // so pfine is in pcoarse and qfine is in qcoarse, pcoarse and qcoarse are different so we create an edge weight
		    if (cg.edge_weights[pcoarse].find(qcoarse)==cg.edge_weights[pcoarse].end())
		      {
			cg.edge_weights[pcoarse][qcoarse] = 0; // zero out weight in new edge
		      }
		      
		    if (pfine==qfine)
		      {
			cg.edge_weights[pcoarse][qcoarse] += vertex_weights[pfine]; // add vertex weight to edge weight
		      }
		    else
		      {
			cg.edge_weights[pcoarse][qcoarse] += edge_weights[pfine][qfine]; // add edge weight to edge weight
		      }
		  }
	      }
      }
    // return result
    return cg;
  }

  void print_info (int verbose=0)
  {
    if (verbose==0) return;
    std::cout << "WeightedGraph [" << std::endl;
    if (verbose>=1)
      std::cout << "graph with " << graph.size() << " vertices" << std::endl;
    if (verbose>=2)
      for (size_t p=0; p<graph.size(); ++p)
	{
	  std::cout << p << "|" << vertex_weights[p] << ": ";
	  for (auto q : graph[p])
	    {
	      std::cout << " " << q << "(" ;
	      if (q!=p) std::cout << edge_weights[p][q];
	      std::cout << ")";
	    }
	  std::cout << std::endl;
	}
    std::cout << "]" << std::endl;
  }
};


/** \brief Describes a decomposition of a grid into subdomains, volume snippets and skeleton snippets
 *
 * \author Peter Bastian
 *
 *
 */
class DomainDecomposition
{
public:
  // the subdomains
  size_t subdomains;                                           // the number of subdomains
  std::vector<size_t> partition;                               // save original partitioning information on finest level for visualization
  std::vector<size_t> vizpartition;                                // partitioning with reduced number of colors
  std::vector<std::set<size_t>> overlappingsubdomains;         // assigns a set of subdomain numbers to each element; mainly used for visualization
  std::shared_ptr<WeightedGraph> p_subdomain_graph;            // weighted graph associated with the domain decomposition on this level using element information

  // volume snippet information
  std::map<std::set<size_t>,size_t> volumesnippet_to_number; // a set of all the snippets we have
  std::vector<std::set<size_t>> number_to_volumesnippet;       // map back volume snippet number to volume snippet
  std::vector<size_t> element_to_volumesnippetnumber;          // assigns number of volume snippet to each mesh element; mainly used for visualization

  // skeleton snippet information
  std::map<std::set<size_t>,size_t> skeletonsnippet_to_number; // a set of all the skeleton snippets
  std::vector<std::set<size_t>> number_to_skeletonsnippet;     // maps back number to skeleton snippet

  // connection of subdomains and snippets on this level
  std::vector<std::vector<size_t>> subdomain_to_volumesnippetnumber; // for each subdomain collect a set of numbers of volume snippets making up that subdomain
  std::vector<std::vector<size_t>> subdomain_to_interior_skeletonsnippetnumber; // for each subdomain store the numbers of skeleton snippets in its interior
  std::vector<std::vector<size_t>> subdomain_to_boundary_skeletonsnippetnumber; // for each subdomain store the numbers of skeleton snippets on its boundary

  // hierarchical information, connection with coarser decomposition; is only filled when this is constructed from a fine domain decomposition
  std::shared_ptr<DomainDecomposition> fine_pdd; // holds pointer to fine domain decomposition if this is not the finest level
  std::vector<std::set<size_t>> subdomain_fine_to_coarse_map;  // subdomain_fine_to_coarse_map[p] contains set of subdomain numbers on this level where p (from fine level) is contained
  std::vector<std::set<size_t>> subdomain_coarse_to_fine_map; // maps subdomain number on the this level to set of subdomain numbers on the finer level that make up the the coarse subdomain
  std::vector<std::set<size_t>> volumesnippet_coarse_to_fine; // volumesnippet_coarse_to_fine[v] gives set of fine volume snippets making up v. Each fine vs contributes to exactly one coarse vs
  std::vector<std::set<size_t>> volumesnippet_coarse_to_fine_skeletonsnippetnumber; // fine skeleton snippets contained in coarse volume snippet
  std::vector<std::set<size_t>> skeletonsnippet_coarse_to_fine; // skeletonsnippet_coarse_to_fine[s] gives set of fine skeleton snippets mapped to s
  
private:

  // maps volume snippet number from fine domain decomposition to this level
  // needs subdomain_fine_to_coarse_map to be filled before
  std::set<size_t> reduce_snippet (const std::set<size_t>& vs)
  {
    std::set<size_t> s;
    for (auto p : vs)
      for (auto q : subdomain_fine_to_coarse_map[p])
	s.insert(q);
    return s;
  }

public:

  template<typename T>
  std::set<T> set_union (const std::set<T>& A, const std::set<T>& B)
  {
    std::set<T> C(A);
    for (const auto& x : B) C.insert(x);
    return C;
  }

  template<typename T>
  std::set<T> set_intersection (const std::set<T>& A, const std::set<T>& B)
  {
    std::set<T> C;
    for (const auto& x : A)
      if (B.count(x)>0)
	C.insert(x);
    return C;
  }

  // construction of empty object
  DomainDecomposition ()
  {}
  
  // create from finer DomainDecomposition and a coarsening of the subdomains
  DomainDecomposition (std::shared_ptr<DomainDecomposition> ptrinputdd, const std::vector<std::set<size_t>>& coarse_to_fine_map)
  {
    Dune::Timer timer;
    timer.reset();
    
    fine_pdd = ptrinputdd; // store pointer to fine decomposition
    DomainDecomposition& inputdd = *ptrinputdd; // a reference to be used within this constructor for ease of writing
    auto nf = inputdd.subdomains;
    auto nc = coarse_to_fine_map.size(); // number of vertices in reduced graph
    subdomains = nc; // this is the number of subdomains

    // we store the coarse to fine map for subdomains
    subdomain_coarse_to_fine_map = coarse_to_fine_map;

    // ... and rebuild fine to coarse map (we had this before ...)
    subdomain_fine_to_coarse_map.resize(nf);
    for (size_t p=0; p<nc; p++) // loop over coarse vertices
      for (auto q : coarse_to_fine_map[p]) // loop over all fine vertices contributing to this coarse vertex
	subdomain_fine_to_coarse_map[q].insert(p);

    // make partition info for visualization
    partition.resize(inputdd.partition.size());
    for (size_t i=0; i<partition.size(); i++)
      for (auto q : subdomain_fine_to_coarse_map[inputdd.partition[i]])
    	{
    	  partition[i] = q; // take the first subdomain
    	  break;
    	}
    std::cout << "XXX == 01 " << timer.elapsed() << std::endl;

    // reduce overlapping subdomain information for the elements
    overlappingsubdomains.resize(inputdd.overlappingsubdomains.size()); // same size because it is for same elements
    // old version:
    // for (size_t i=0; i<inputdd.overlappingsubdomains.size(); i++)
    //   overlappingsubdomains[i] = reduce_snippet(inputdd.overlappingsubdomains[i]);
    // new version:
    std::vector<std::set<size_t>> finevsnumber_to_coarsevs(inputdd.volumesnippet_to_number.size());
    for (auto it=inputdd.volumesnippet_to_number.begin(); it!=inputdd.volumesnippet_to_number.end(); ++it) // iterate over all volume snippets in input
      finevsnumber_to_coarsevs[it->second] = reduce_snippet(it->first); // map snippet from fine to coarse and create entry
    for (size_t i=0; i<inputdd.overlappingsubdomains.size(); i++)
      overlappingsubdomains[i] = finevsnumber_to_coarsevs[inputdd.element_to_volumesnippetnumber[i]];
    std::cout << "XXX == 02 " << timer.elapsed() << std::endl;

    // construct subdomain graph from element information
    p_subdomain_graph = std::shared_ptr<WeightedGraph>(new WeightedGraph()); // allocate graph
    WeightedGraph& g=*p_subdomain_graph; // for simplicity
    g.graph.resize(subdomains);
    g.vertex_weights.resize(subdomains,0);
    g.edge_weights.resize(subdomains);
    for (auto& S : overlappingsubdomains) // loop over subdomain information of all elements
      for (auto p : S)
	{
	  g.vertex_weights[p]++; // increment vertex weight by one element
	  for (auto q : S)
	    {
	      g.graph[p].insert(q); // couple the two subdomains because they share an element
	      if (p!=q) // edge weights are only calculated for offdiagonal edges
		if (g.edge_weights[p].find(q)==g.edge_weights[p].end())
		  g.edge_weights[p][q] = 1; // new entry
		else
		  g.edge_weights[p][q] += 1; // increment existing entry
	    }
	}
    // analyse subdomain graph
    size_t min_edge_weight = overlappingsubdomains.size();
    for (size_t p=0; p<g.edge_weights.size(); p++) // loop over subdomains
      for (auto& entry : g.edge_weights[p]) // loop over all neighbors of p
	min_edge_weight = std::min(entry.second,min_edge_weight);
    // for (size_t p=0; p<g.edge_weights.size(); p++) // loop over subdomains
    //   for (auto& entry : g.edge_weights[p]) // loop over all neighbors of p
    // 	if (entry.second<=2*min_edge_weight && entry.first>p)
    // 	  std::cout << "subdomain " << p << " and " << entry.first << " have " << entry.second << " joint elements" << std::endl;
    std::cout << "XXX == 03 " << timer.elapsed() << std::endl;

    // color map for the partitioning
    vizpartition.resize(partition.size()); // to be exported
    std::vector<size_t> color(subdomains); // map subdomains to colors
    for (size_t p=0; p<subdomains; p++)
      {
	std::set<size_t> usedcolors;
	for (auto q : p_subdomain_graph->graph[p])
	  if (q<p)
	    usedcolors.insert(color[q]);
	for (size_t c=0; c<10000; c++)
	  if (usedcolors.count(c)==0)
	    {
	      color[p] = c;
	      break;
	    }
      }
    for (size_t i=0; i<partition.size(); i++)
      vizpartition[i] = color[partition[i]];
    std::cout << "XXX == 04 " << timer.elapsed() << std::endl;

    // collect the new volume snippets
    for (auto it=inputdd.volumesnippet_to_number.begin(); it!=inputdd.volumesnippet_to_number.end(); ++it) // iterate over all volume snippets in input
      volumesnippet_to_number[reduce_snippet(it->first)] = 0; // map snippet from fine to coarse and create entry
    // now assign numbering of snippets
    for (auto it=volumesnippet_to_number.begin(); it!=volumesnippet_to_number.end(); ++it)
      {
        number_to_volumesnippet.push_back(it->first);
        it->second = number_to_volumesnippet.size()-1;
      }
    std::cout << "XXX == 05 " << timer.elapsed() << std::endl;

    // build hierarchy of volume snippets
    volumesnippet_coarse_to_fine.resize(number_to_volumesnippet.size()); // resize to number of volume snippets on this level
    for (auto it=inputdd.volumesnippet_to_number.begin(); it!=inputdd.volumesnippet_to_number.end(); ++it) // iterate over all volume snippets in input
      volumesnippet_coarse_to_fine[volumesnippet_to_number[reduce_snippet(it->first)]].insert(it->second);
    std::cout << "XXX == 06 " << timer.elapsed() << std::endl;

    // store snippet number for each element
    // this is to save memory and time & we may double check the computation of the volume snippets
    std::vector<size_t> finevsnumber_to_coarsevsnumber(inputdd.volumesnippet_to_number.size());
    for (size_t i=0; i<finevsnumber_to_coarsevs.size(); i++) // iterate over all volume snippets in input
      finevsnumber_to_coarsevsnumber[i] = volumesnippet_to_number[finevsnumber_to_coarsevs[i]];
    element_to_volumesnippetnumber.resize(overlappingsubdomains.size());
    for (size_t i=0; i<element_to_volumesnippetnumber.size(); i++)
      element_to_volumesnippetnumber[i] = finevsnumber_to_coarsevsnumber[inputdd.element_to_volumesnippetnumber[i]];
      // {
      // 	auto it = volumesnippet_to_number.find(overlappingsubdomains[i]);
      // 	if (it!=volumesnippet_to_number.end())  
      // 	  element_to_volumesnippetnumber[i] = it->second;
      // 	else
      // 	  DUNE_THROW(Dune::Exception, "volumesnippet_to_number access failed");
      // }
    std::cout << "XXX == 07 " << timer.elapsed() << std::endl;

    // now for the skeleton snippets
    for (auto& skeletonsnippet : inputdd.number_to_skeletonsnippet)
      {
        auto it = skeletonsnippet.begin(); // points to first snippet number
        size_t left_volumesnippetnumber = *it;
        ++it;
        size_t right_volumesnippetnumber = *it; // so now we have the numbers of the two volume snippets to left and right
        auto& left_volumesnippet = inputdd.number_to_volumesnippet[left_volumesnippetnumber];
	auto left_reduced_snippet = reduce_snippet(left_volumesnippet);
        auto& right_volumesnippet = inputdd.number_to_volumesnippet[right_volumesnippetnumber];
	auto right_reduced_snippet = reduce_snippet(right_volumesnippet);
	if (left_reduced_snippet!=right_reduced_snippet)
	  {
	    // we have a snippet
	    std::set<size_t> key; // a key to index the interface between two snippets
	    key.insert(volumesnippet_to_number[left_reduced_snippet]);
	    key.insert(volumesnippet_to_number[right_reduced_snippet]);
	    skeletonsnippet_to_number[key] = 0;
	  }
      }
    for (auto it=skeletonsnippet_to_number.begin(); it!=skeletonsnippet_to_number.end(); ++it)
      {
        number_to_skeletonsnippet.push_back(it->first);
        it->second = number_to_skeletonsnippet.size()-1;
      }
    std::cout << "XXX == 08 " << timer.elapsed() << std::endl;

    // build hierarchy for skeleton snippets
    volumesnippet_coarse_to_fine_skeletonsnippetnumber.resize(number_to_volumesnippet.size()); // resize to number of volume snippets on this level
    skeletonsnippet_coarse_to_fine.resize(number_to_skeletonsnippet.size()); // resize to number of skeleton snippets on this level
    for (size_t i=0;  i<inputdd.number_to_skeletonsnippet.size(); ++i) // loop over all fine skeleton snippets
      {
	auto& skeletonsnippet = inputdd.number_to_skeletonsnippet[i]; // fine skeleton snippet with number i
        auto it = skeletonsnippet.begin(); // points to first v snippet number
        size_t left_volumesnippetnumber = *it;
        ++it;
        size_t right_volumesnippetnumber = *it; // so now we have the numbers of the two volume snippets to left and right
	auto left_reduced_volumesnippet_number = volumesnippet_to_number[reduce_snippet(inputdd.number_to_volumesnippet[left_volumesnippetnumber])];
	auto right_reduced_volumesnippet_number = volumesnippet_to_number[reduce_snippet(inputdd.number_to_volumesnippet[right_volumesnippetnumber])];
	if (left_reduced_volumesnippet_number!=right_reduced_volumesnippet_number)
	  {
	    // to left and right are different volume snippets on the coarse level
	    std::set<size_t> key; // a key to index the interface between two snippets
	    key.insert(left_reduced_volumesnippet_number);
	    key.insert(right_reduced_volumesnippet_number);
	    skeletonsnippet_coarse_to_fine[skeletonsnippet_to_number[key]].insert(i);
	  }
	else
	  {
	    // to left and right are the same volume snippets on the coarse level
	    // so this fine level skeleton snippet is inside a coarse volume snippet
	    volumesnippet_coarse_to_fine_skeletonsnippetnumber[left_reduced_volumesnippet_number].insert(i);
	  }
      }
    std::cout << "XXX == 09 " << timer.elapsed() << std::endl;

    // subdomains and snippets;
    // determine volume snippets making up each subdomain
    subdomain_to_volumesnippetnumber.resize(subdomains); // for each subdomain collect a set of snippet numbers making up that subdomain
    for (const auto& entry : volumesnippet_to_number) // loop over all volume snippets
      for (auto p : entry.first) // loop over all subdomains the volume snippet is in
	subdomain_to_volumesnippetnumber[p].push_back(entry.second); // so subdomain p contains that snippet
    std::cout << "XXX == 10 " << timer.elapsed() << std::endl;

    // skeleton snippets contributing to each subdomains
    subdomain_to_interior_skeletonsnippetnumber.resize(subdomains); // for each subdomain store the numbers of skeleton snippets in its interior
    subdomain_to_boundary_skeletonsnippetnumber.resize(subdomains); // for each subdomain store the numbers of skeleton snippets on its boundary
    for (const auto& entry : skeletonsnippet_to_number) // loop over all skeleton snippets
      {
        auto it = entry.first.begin(); // points to first snippet number
        auto left_volumesnippetnumber = *it;
        auto& left_volumesnippet = number_to_volumesnippet[left_volumesnippetnumber];
        ++it;
        auto right_volumesnippetnumber = *it; // so now we have the numbers of the two volume snippets to left and right
        auto& right_volumesnippet = number_to_volumesnippet[right_volumesnippetnumber];
        for (auto p : left_volumesnippet)
          if (right_volumesnippet.count(p)>0)
            subdomain_to_interior_skeletonsnippetnumber[p].push_back(entry.second); // p contains left and right snippet so it is interior
          else
            subdomain_to_boundary_skeletonsnippetnumber[p].push_back(entry.second); // p contains left snippet but not right
        for (auto q : right_volumesnippet)
          if (left_volumesnippet.count(q)==0)
            subdomain_to_boundary_skeletonsnippetnumber[q].push_back(entry.second); // q contains right snippet but not left
      }
    std::cout << "XXX == 11 " << timer.elapsed() << std::endl;
  }
  
  // construct from a grid level (typically the finest level)
  template<typename GV>
  DomainDecomposition (const GV& gv, size_t subdomains_, std::vector<size_t> parts, size_t overlap, std::string extensionmethod, size_t drop_small_overlap, bool verbose=false)
  {
    // store #subdomains
    subdomains = 1; for (auto i : parts) subdomains *= i;
    if (subdomains!=subdomains_)
      DUNE_THROW(Dune::Exception, "coordinate_partitioning: number of subdomains in parts does not match number of subdomains");

    std::cout << "before partitioning into " << subdomains << " subdomains" << std::endl;
    
    // partition the elements using ParMetis; assigns one number to each element
    partition = coordinate_partitioning(gv,parts);

    std::cout << "I have partitioned" << std::endl;

    build0(gv,overlap,extensionmethod,drop_small_overlap,verbose);
  }

  // construct from a grid level (typically the finest level)
  template<typename GV>
  DomainDecomposition (const GV& gv, MPI_Comm comm, size_t subdomains_, size_t overlap, std::string extensionmethod, size_t drop_small_overlap, bool verbose=false)
    : subdomains(subdomains_)
  {
    // partition the elements using ParMetis; assigns one number to each element
    partition = parmetis_partitioning(gv,comm,subdomains);

    build0(gv,overlap,extensionmethod,drop_small_overlap,verbose);
  }

private:
  // construct from a grid level (typically the finest level)
  template<typename GV>
  void build0 (const GV& gv, size_t overlap, std::string extensionmethod, size_t drop_small_overlap, bool verbose)
  {

    // generate overlap: need to decide on size of overlap as well as extension method; assigns a set of numbers to each element
    Dune::Timer timer;
    timer.reset();
    overlappingsubdomains = grow_subdomains(gv,subdomains,partition,overlap,extensionmethod);
    std::cout << "XXX DD 02 " << timer.elapsed() << std::endl;

    // grid information
    const int dim = GV::dimension;
    auto& indexset = gv.indexSet(); // to attach data to elements

    // construct subdomain graph from element information
    p_subdomain_graph = std::shared_ptr<WeightedGraph>(new WeightedGraph()); // allocate graph
    WeightedGraph& g=*p_subdomain_graph; // for simplicity
    g.graph.resize(subdomains);
    g.vertex_weights.resize(subdomains,0);
    g.edge_weights.resize(subdomains);
    for (auto& S : overlappingsubdomains) // loop over subdomain information of all elements
      for (auto p : S)
	{
	  g.vertex_weights[p]++; // increment vertex weight by one element
	  for (auto q : S)
	    {
	      g.graph[p].insert(q); // couple the two subdomains because they share an element
	      if (p!=q) // edge weights are only calculated for offdiagonal edges
		if (g.edge_weights[p].find(q)==g.edge_weights[p].end())
		  g.edge_weights[p][q] = 1; // new entry
		else
		  g.edge_weights[p][q] += 1; // increment existing entry
	    }
	}
    // analyse subdomain graph
    size_t min_edge_weight = overlappingsubdomains.size();
    for (size_t p=0; p<g.edge_weights.size(); p++) // loop over subdomains
      for (auto& entry : g.edge_weights[p]) // loop over all neighbors of p
	min_edge_weight = std::min(entry.second,min_edge_weight);
    // for (size_t p=0; p<g.edge_weights.size(); p++) // loop over subdomains
    //   for (auto& entry : g.edge_weights[p]) // loop over all neighbors of p
    // 	if (entry.second<=10 && entry.first>p)
    // 	  std::cout << "subdomain " << p << " and " << entry.first << " have " << entry.second << " joint elements" << std::endl;
    std::cout << "XXX DD 03 " << timer.elapsed() << std::endl;

    // try to repair the domain decomposition
    for (size_t p=0; p<g.edge_weights.size(); p++) // loop over subdomains
      for (auto& entry : g.edge_weights[p]) // loop over all neighbors of p
	if (p!=entry.first && entry.second<drop_small_overlap)
	  {
	    // now we have at most 4 joint elements between p and q
	    auto q = entry.first;

	    // lets try to find them one by one and remove p and q
	    for (size_t i=0; i<overlappingsubdomains.size(); i++) // loop over subdomain information of all elements
	      {
		auto& S = overlappingsubdomains[i];
		if (S.count(p)>0 && S.count(q)>0)
		  if (S.size()>2)
		    {
		      // now there is at least one other subdomain for this element, lets remove p and q
		      //std::cout << vs_to_str(S) << " -> ";
		      S.erase(p);
		      S.erase(q);
		      //std::cout << vs_to_str(S) << std::endl;
		    }
		  else
		    {
		      std::cout << vs_to_str(S) << " cannot remove " << p << " and " << q << std::endl;
		      DUNE_THROW(Dune::Exception, "cannot repair domain decomposition");
		    }
	      }
	  }
    std::cout << "XXX DD 04 " << timer.elapsed() << std::endl;
    
    // color map for the partitioning
    vizpartition.resize(partition.size()); // to be exported
    std::vector<size_t> color(subdomains); // map subdomains to colors
    for (size_t p=0; p<subdomains; p++)
      {
	std::set<size_t> usedcolors;
	for (auto q : p_subdomain_graph->graph[p])
	  if (q<p)
	    usedcolors.insert(color[q]);
	for (size_t c=0; c<10000; c++)
	  if (usedcolors.count(c)==0)
	    {
	      color[p] = c;
	      break;
	    }
      }
    for (size_t i=0; i<partition.size(); i++)
      vizpartition[i] = color[partition[i]];
    std::cout << "XXX DD 05 " << timer.elapsed() << std::endl;
    
    // generate volume snippets and subdomains as unions of snippets
    // a volume snippet is a set of elements together with local<->global numbering and stiffness matrix
    // a volume snippet is represented by its key: an ordered list of numbers of subdomains and each snippet gets a number
    for (auto& S : overlappingsubdomains)
      volumesnippet_to_number[S] = 0; // do not know the number yet
    if (verbose) std::cout << volumesnippet_to_number.size() << " volume snippets found" << std::endl;
    for (auto it=volumesnippet_to_number.begin(); it!=volumesnippet_to_number.end(); ++it)
      {
        number_to_volumesnippet.push_back(it->first);
        it->second = number_to_volumesnippet.size()-1;
      }

    // store snippet number for each element
    // this is to save memory and time
    element_to_volumesnippetnumber.resize(overlappingsubdomains.size());
    for (size_t i=0; i<element_to_volumesnippetnumber.size(); i++)
      element_to_volumesnippetnumber[i] = volumesnippet_to_number[overlappingsubdomains[i]];

    // determine snippets making up each subdomain
    subdomain_to_volumesnippetnumber.resize(subdomains); // for each subdomain collect a set of snippet numbers making up that subdomain
    for (const auto& entry : volumesnippet_to_number) // loop over all volume snippets
      for (auto p : entry.first) // loop over all subdomains the volume snippet is in
	subdomain_to_volumesnippetnumber[p].push_back(entry.second); // so subdomain p contains that snippet
    // and make some statistics
    size_t total_snippets = 0;
    for (size_t p=0; p<subdomains; p++)
      total_snippets += subdomain_to_volumesnippetnumber[p].size();
    if (verbose) std::cout << ((double)total_snippets)/subdomains << " volume snippets per subdomain" << std::endl;

    // now for the skeleton snippets
    // a skeleton snippet is a set of faces seperating two volume snippets
    // a skeleton snippet is represented by a pair of volume snippet numbers
    for (const auto& e : elements(gv)) // loop over elements in the mesh
      for (const auto& is : intersections(gv,e)) // loop over faces
        if (is.neighbor()) // if this is an intersection with two elements
          if (element_to_volumesnippetnumber[indexset.index(e)]!=element_to_volumesnippetnumber[indexset.index(is.outside())]) // if this face seperates two different snippets
            {
              // then we have to store the data for this intersection seperately
              std::set<size_t> key; // a key to index the interface between two snippets
              key.insert(element_to_volumesnippetnumber[indexset.index(e)]);
              key.insert(element_to_volumesnippetnumber[indexset.index(is.outside())]);
              skeletonsnippet_to_number[key] = 0;
            }
    if (verbose) std::cout << skeletonsnippet_to_number.size() << " skeleton snippets found" << std::endl;
    for (auto it=skeletonsnippet_to_number.begin(); it!=skeletonsnippet_to_number.end(); ++it)
      {
        number_to_skeletonsnippet.push_back(it->first);
        it->second = number_to_skeletonsnippet.size()-1;
      }

    // subdomains to skeleton snippets
    subdomain_to_interior_skeletonsnippetnumber.resize(subdomains); // for each subdomain store the numbers of skeleton snippets in its interior
    subdomain_to_boundary_skeletonsnippetnumber.resize(subdomains); // for each subdomain store the numbers of skeleton snippets on its boundary
    for (const auto& entry : skeletonsnippet_to_number) // loop over all skeleton snippets
      {
        auto it = entry.first.begin(); // points to first snippet number
        auto left_volumesnippetnumber = *it;
        auto& left_volumesnippet = number_to_volumesnippet[left_volumesnippetnumber];
        ++it;
        auto right_volumesnippetnumber = *it; // so now we have the numbers of the two volume snippets to left and right
        auto& right_volumesnippet = number_to_volumesnippet[right_volumesnippetnumber];
        for (auto p : left_volumesnippet)
          if (right_volumesnippet.count(p)>0)
            subdomain_to_interior_skeletonsnippetnumber[p].push_back(entry.second); // p contains left and right snippet so it is interior
          else
            subdomain_to_boundary_skeletonsnippetnumber[p].push_back(entry.second); // p contains left snippet but not right
        for (auto q : right_volumesnippet)
          if (left_volumesnippet.count(q)==0)
            subdomain_to_boundary_skeletonsnippetnumber[q].push_back(entry.second); // q contains right snippet but not left
      }
    std::cout << "XXX DD 06 " << timer.elapsed() << std::endl;
  }

public:
  std::string vs_to_str (std::set<size_t> vs)
  {
    std::ostringstream s;
    s << "(";
    for (auto i : vs) s << " " << i;
    s << " )";
    return s.str();
  }

  std::string ss_to_str (std::set<size_t> ss)
  {
    auto it = ss.begin(); // points to first snippet number
    auto left_volumesnippetnumber = *it;
    auto& left_volumesnippet = number_to_volumesnippet[left_volumesnippetnumber];
    ++it;
    auto right_volumesnippetnumber = *it; // so now we have the numbers of the two volume snippets to left and right
    auto& right_volumesnippet = number_to_volumesnippet[right_volumesnippetnumber];
    std::ostringstream s;
    s << "[";
    s << vs_to_str(left_volumesnippet) << "|" << vs_to_str(right_volumesnippet);
    s << " ]";
    return s.str();
  }

  void print_info (int verbose=0)
  {
    // list
    if (verbose==0)
      return;

    std::cout << "DomainDecomposition [" << std::endl;
    if (verbose>=1)
      {
	std::cout << number_to_volumesnippet.size() << " volume snippets" << std::endl;
	std::cout << number_to_skeletonsnippet.size() << " skeleton snippets" << std::endl;
	std::cout << subdomains << " subdomains" << std::endl;
      }
    if (verbose>=2)
      {
	for (size_t p = 0; p<subdomains; ++p)
	  {
	    std::cout << p << ": " << std::endl;
	    std::cout << " volume snippets:";
	    for (auto i : subdomain_to_volumesnippetnumber[p])
	      std::cout << " " << vs_to_str(number_to_volumesnippet[i]);
	    std::cout << std::endl;
	    std::cout << " interior skeleton snippets:";
	    for (auto i : subdomain_to_interior_skeletonsnippetnumber[p])
	      std::cout << " " << ss_to_str(number_to_skeletonsnippet[i]);
	    std::cout << std::endl;
	    std::cout << " boundary skeleton snippets:";
	    for (auto i : subdomain_to_boundary_skeletonsnippetnumber[p])
	      std::cout << " " << ss_to_str(number_to_skeletonsnippet[i]);
	    std::cout << std::endl;
	  }
	for (size_t i=0; i<volumesnippet_coarse_to_fine.size(); ++i)
	  {
	    std::cout << "fine vs contributing to vs " << i  << " " << vs_to_str(number_to_volumesnippet[i]) << ": ";
	    for (auto j : volumesnippet_coarse_to_fine[i]) std::cout << " " << fine_pdd->vs_to_str(fine_pdd->number_to_volumesnippet[j]);
	    std::cout << std::endl;
	    std::cout << "fine ss interior to vs " << i  << " " << vs_to_str(number_to_volumesnippet[i]) << ": ";
	    for (auto j : volumesnippet_coarse_to_fine_skeletonsnippetnumber[i]) std::cout << " " << fine_pdd->ss_to_str(fine_pdd->number_to_skeletonsnippet[j]);
	    std::cout << std::endl;
	  }
	for (size_t i=0; i<skeletonsnippet_coarse_to_fine.size(); ++i)
	  {
	    std::cout << "fine ss contributing to ss " << i  << " " << ss_to_str(number_to_skeletonsnippet[i]) << ": ";
	    for (auto j : skeletonsnippet_coarse_to_fine[i]) std::cout << " " << fine_pdd->ss_to_str(fine_pdd->number_to_skeletonsnippet[j]);
	    std::cout << std::endl;							       
	  }
      }
    std::cout << "]" << std::endl;
  }
}; // class




/** \brief manage degrees of freedom to subdomains and snippets
 *
 * \author Peter Bastian
 *
 *
 */
class DomainDecompositionDOFMapper
{
public:
  size_t n_global_dofs;  // the global number of dof (bolcks) on this level; use this to allocate ISTL vector
  size_t blocksize;      // size of blocks in ISTL vector; always 1 one coarse levels
  
  std::shared_ptr<DomainDecomposition> pdd;                                 // associated domain decomposition
  std::shared_ptr<DomainDecompositionDOFMapper> fine_pdddm;                 // DomainDecompositionDOFMapper on next finer level if this is a coarse level
  
  std::vector<std::vector<size_t>> volumesnippetnumber_local_to_global;     // volumesnippetnumber_local_to_global[s][i] is global dof (block) number of dof (block) i in volume snippet number s
  std::vector<std::map<size_t,size_t>> volumesnippetnumber_global_to_local; // volumesnippetnumber_global_to_local[i][s] gives the local index in volume snippet s for global index i 

  std::vector<std::vector<size_t>> skeletonsnippetnumber_local_to_global;     // skeletonsnippetnumber_local_to_global[s][i] is global dof (block) number of dof (block) i in skeleton snippet number s
  std::vector<std::map<size_t,size_t>> skeletonsnippetnumber_global_to_local; // skeletonsnippetnumber_global_to_local[i][s] gives the local index in skeleton snippet s for global index i 

  std::vector<std::vector<size_t>> subdomainnumber_local_to_global;     // subdomainnumber_local_to_global[s][i] is global dof (block) number of dof (block) i in subdomain s
  std::vector<std::map<size_t,size_t>> subdomainnumber_global_to_local; // subdomainnumber_global_to_local[i][s] gives the local index in subdomain s for global index i
  std::vector<std::vector<bool>> subdomain_boundary; // boundary[p][i] is true if dof (block) i in subdomain p is at the subdomain boundary ... only valid on coarse levels so far

  std::shared_ptr<WeightedGraph> p_subdomain_graph; // weighted graph associated with the domain decomposition on this level using dof information

  // local_to_global mapping for basis vectors from the subdomains on this level;
  // this determines to global dofs on the coarse grid
  std::vector<size_t> n_basis_vectors_per_subdomain; // information given in from outside
  std::vector<bool> disconnected_subdomains; // disconnected_subdomains[p] true when this subdomain is not connected (has >1 zero EV)
  std::vector<std::vector<size_t>> basisvector_local_to_global; // fine_subdomainnumber_local_to_global[p][k] gives global number of basisvector k in subdomain p ON FINE LEVEL
  std::vector<size_t> basisvector_global_to_subdomainnumber;  //there is a 1 to 1 correspondance between basis vectors in sudomains on fine level with global dofs on this level
  std::vector<size_t> basisvector_global_to_local;

  // Constructor for finest level: dofs are determined by a grid function space
  template<typename GFS>
  DomainDecompositionDOFMapper (const GFS& gfs, std::shared_ptr<DomainDecomposition> pdd_)
    : pdd(pdd_)
  {
    // some types
    using GV = typename GFS::Traits::GridView;
    using IS = typename GV::IndexSet;
    using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
    using LFSCache = Dune::PDELab::LFSIndexCache<LFS>;
    using RF = size_t;
    using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
    using ISTLV = Dune::PDELab::Backend::Native<Z>;
    using VIEW = typename Z::template LocalView<LFSCache>;

    // local variables
    const GV& gv(gfs.gridView()); // from the gfs
    const IS& indexset(gv.indexSet()); // from the gv
    LFS lfs(gfs); // a local function space for the gfs
    LFSCache lfscache(lfs); // lfs cache for accessing dofs of an element
    Z pdelab_dofnumber(gfs); // one-dimensional flat index for each dof
    ISTLV& dofnumber = Dune::PDELab::Backend::native(pdelab_dofnumber); // a vector storing the global dof number in each dof
    VIEW view(pdelab_dofnumber); // local part of the vector associated with an element
    auto& dd = *pdd; // use reference, its nicer

    // first we fill a vector with the global dof numbers
    // we assume that our ISTL vector is a BlockVector<FieldVector<>>
    RF dof_counter = 0; 
    for (int i=0; i<dofnumber.size(); i++)
      for (int j=0; j<dofnumber[i].size(); j++)
        dofnumber[i][j] = dof_counter++;

    // determine the block size
    blocksize = dofnumber[0].size();

    // now we can associate global dofs with volume snippets;  collect global dofs for each skeleton snippet first
    size_t n_volume_snippets = dd.number_to_volumesnippet.size();
    {
      std::vector<std::set<size_t>> volumesnippetnumber_to_globaldofs(n_volume_snippets); // volumesnippetnumber_to_globaldofs[i] gives set of dof numbers for volume snippet number i
      for (const auto& e : elements(gv)) // loop over all elements
	{
	  // read all dof numbers associated with this element
	  lfs.bind(e);
	  lfscache.update();
	  std::vector<RF> localdata(lfs.size());
	  view.bind(lfscache);
	  view.read(localdata);
	  view.unbind();
	  
	  // scatter dofs to volume snippets
	  auto number = dd.element_to_volumesnippetnumber[indexset.index(e)]; // the volume snippet number for this element
	  for (size_t i = 0; i<lfs.size(); ++i) // all these global dofs are in that volume snippet
	    volumesnippetnumber_to_globaldofs[number].insert(localdata[i]/blocksize);
	}
      
      // now we can make local to global map for each volume snippet
      volumesnippetnumber_local_to_global.resize(n_volume_snippets); // resize for fill-in
      for (size_t i=0; i<n_volume_snippets; i++)
	for (auto g : volumesnippetnumber_to_globaldofs[i])
	  volumesnippetnumber_local_to_global[i].push_back(g); // add global dof number to volume snippet number i
      
      // invert the local to global to get global to local
      n_global_dofs = dofnumber.N(); // number of global dof blocks
      volumesnippetnumber_global_to_local.resize(n_global_dofs);
      for (size_t s=0; s<n_volume_snippets; s++)
	for (size_t i=0; i<volumesnippetnumber_local_to_global[s].size(); i++)
	  volumesnippetnumber_global_to_local[volumesnippetnumber_local_to_global[s][i]][s] = i;
    }
    
    // now for the skeleton snippets; collect global dofs for each skeleton snippet first
    size_t n_skeleton_snippets = dd.number_to_skeletonsnippet.size();
    {
      std::vector<std::set<size_t>> skeletonsnippetnumber_to_globaldofs(n_skeleton_snippets); // skeletonsnippetnumber_to_globaldofs[i] gives set of dof numbers for skeleton snippet number i
      for (const auto& e : elements(gv)) // loop over elements in the mesh
	{
	  // read all dof numbers associated with this element
	  lfs.bind(e);
	  lfscache.update();
	  std::vector<RF> localdata_inside(lfs.size());
	  view.bind(lfscache);
	  view.read(localdata_inside);
	  view.unbind();

	  for (const auto& is : intersections(gv,e)) // loop over faces
	    if (is.neighbor()) // if this is an interior intersection with two elements
	      if (dd.element_to_volumesnippetnumber[indexset.index(e)]!=dd.element_to_volumesnippetnumber[indexset.index(is.outside())]) // if this face seperates two different snippets
		{
		  // now this intersection is on the interface beween two snippets and we have to assemble its contribution seperately
		  std::set<size_t> key; // a key to index the interface between two snippets
		  key.insert(dd.element_to_volumesnippetnumber[indexset.index(e)]);
		  key.insert(dd.element_to_volumesnippetnumber[indexset.index(is.outside())]);
		  auto number = dd.skeletonsnippet_to_number[key];

		  // read all dof numbers associated with outside element
		  lfs.bind(is.outside());
		  lfscache.update();
		  std::vector<RF> localdata_outside(lfs.size());
		  view.bind(lfscache);
		  view.read(localdata_outside);
		  view.unbind();

		  // now put all dofs of both elements in skeleton snippet
		  for (size_t i = 0; i<lfs.size(); ++i)
		    {
		      skeletonsnippetnumber_to_globaldofs[number].insert(localdata_inside[i]/blocksize);
		      skeletonsnippetnumber_to_globaldofs[number].insert(localdata_outside[i]/blocksize);
		    }
		}
	}

      // now we can make local to global map for each skeleton snippet
      skeletonsnippetnumber_local_to_global.resize(n_skeleton_snippets); // resize for fill-in
      for (size_t i=0; i<n_skeleton_snippets; i++)
	for (auto g : skeletonsnippetnumber_to_globaldofs[i])
	  skeletonsnippetnumber_local_to_global[i].push_back(g); // add global dof number to volume snippet number i

      // invert the local to global to get global to local
      skeletonsnippetnumber_global_to_local.resize(dofnumber.N()); // n_globa_dofs is still the same ...
      for (size_t s=0; s<n_skeleton_snippets; s++)
	for (size_t i=0; i<skeletonsnippetnumber_local_to_global[s].size(); i++)
	  skeletonsnippetnumber_global_to_local[skeletonsnippetnumber_local_to_global[s][i]][s] = i;
    }

    // now the same thing for the subdomain local <-> global
    size_t n_subdomains = dd.subdomains;
    subdomainnumber_local_to_global.resize(n_subdomains); // resize for fill-in
    subdomainnumber_global_to_local.resize(dofnumber.N()); // n_global_dofs is still the same ...
    for (size_t p=0; p<n_subdomains; ++p) // loop over all subdomains
      {
	std::set<size_t> globaldofs; // find global dofs for this subdomain
	
	for (auto volumesnippetnumber : dd.subdomain_to_volumesnippetnumber[p]) // loop over all volumesnippets in that subdomain
	  for (auto g : volumesnippetnumber_local_to_global[volumesnippetnumber]) // loop over all global dof numbers in that snippet
	    globaldofs.insert(g); // add global dof to subdomain

	// now we can build up local to global map
	for (auto g : globaldofs)
	  subdomainnumber_local_to_global[p].push_back(g);

	// and global to local
	for (size_t i=0; i<subdomainnumber_local_to_global[p].size(); i++)
	  subdomainnumber_global_to_local[subdomainnumber_local_to_global[p][i]][p] = i;
      }    

    // construct weighted subdomain graph on this level
    // weight of a subdomain: number of dofs in this subdomain
    // weight of an edge (p,q): number of dofs shared by p and q
    p_subdomain_graph = std::shared_ptr<WeightedGraph>(new WeightedGraph()); // allocate graph
    WeightedGraph& g=*p_subdomain_graph; // for simplicity
    g.graph.resize(pdd->subdomains);
    g.vertex_weights.resize(pdd->subdomains,0);
    g.edge_weights.resize(pdd->subdomains);
    for (size_t j=0; j<subdomainnumber_global_to_local.size(); j++) // loop over global dof numbers
      for (auto& itp : subdomainnumber_global_to_local[j]) // loop over subdomain numbers sharing this degree of freedom
	{
	  g.vertex_weights[itp.first]++; // increment vertex weight by one dof
	  for (auto& itq : subdomainnumber_global_to_local[j]) // loop over subdomain numbers sharing this degree of freedom
	    {
	      g.graph[itp.first].insert(itq.first); // couple the two subdomains because they share a degree of freedom
	      if (itp.first!=itq.first) // edge weights are only calculated for offdiagonal edges
		if (g.edge_weights[itp.first].find(itq.first)==g.edge_weights[itp.first].end())
		  g.edge_weights[itp.first][itq.first] = 1; // new entry
		else
		  g.edge_weights[itp.first][itq.first] += 1; // increment existing entry
	    }
	}
    // analyse subdomain graph
    size_t min_edge_weight = subdomainnumber_global_to_local.size();
    for (size_t p=0; p<g.edge_weights.size(); p++) // loop over subdomains
      for (auto& entry : g.edge_weights[p]) // loop over all neighbors of p
	min_edge_weight = std::min(entry.second,min_edge_weight);
    // for (size_t p=0; p<g.edge_weights.size(); p++) // loop over subdomains
    //   for (auto& entry : g.edge_weights[p]) // loop over all neighbors of p
    // 	if (entry.second<=2*min_edge_weight && entry.first>p)
    // 	  std::cout << "subdomain " << p << " and " << entry.first << " have " << entry.second << " joint dofs" << std::endl;
  }


  // fill in additional information about basis vectors on a level
  // this can only be done after the eigenproblem has been solved
  void fill_basisvector_information (const std::vector<size_t>& n_basis_vectors_per_subdomain_, const std::vector<bool>& disconnected_subdomains_)
  {
    if (pdd->subdomains!=n_basis_vectors_per_subdomain_.size())
      DUNE_THROW(Dune::Exception, "DomainDecompositionDOFMapper::fill_basisvector_information number of subdomains on fine level mismatch");
    if (pdd->subdomains!=disconnected_subdomains_.size())
      DUNE_THROW(Dune::Exception, "DomainDecompositionDOFMapper::fill_basisvector_information number of subdomains on fine level mismatch");

    // store number of basis vectors computed by eigensolver
    n_basis_vectors_per_subdomain = n_basis_vectors_per_subdomain_;
    disconnected_subdomains = disconnected_subdomains_;

    // basis vector local to global on this level
    basisvector_local_to_global.resize(n_basis_vectors_per_subdomain.size());
    for (size_t p=0; p<n_basis_vectors_per_subdomain.size(); ++p) // loop over all subdomains
      {
	basisvector_local_to_global[p].resize(n_basis_vectors_per_subdomain[p]); // make space	
	for (size_t k=0; k<n_basis_vectors_per_subdomain[p]; ++k)
	  {
	    basisvector_global_to_subdomainnumber.push_back(p);
	    basisvector_global_to_local.push_back(k);
	    basisvector_local_to_global[p][k] = basisvector_global_to_local.size()-1;
	  }
      }
  }

  // Constructor for a coarse level: dofs are determined by basisvectors in subdomains on next finer level
  DomainDecompositionDOFMapper (std::shared_ptr<DomainDecomposition> pdd_, // domain decomposition on *this* level which points also to next finer level
				std::shared_ptr<DomainDecompositionDOFMapper> fine_pdddm_ // dof information on next finer level
			       )
    : pdd(pdd_), fine_pdddm(fine_pdddm_)
  {
    // general information
    blocksize = 1; // always 1 on coarse level
    n_global_dofs = fine_pdddm->basisvector_global_to_local.size();

    // now associate global dofs on this level with volume snippets on this level
    // this is independent of the domain decomposition on this level but only depends on domain decomposition of fine level
    // approach: for all vs on this level loop over vs on fine level contributing to it and determine the subdomains it is in
    volumesnippetnumber_local_to_global.resize(pdd->number_to_volumesnippet.size()); // how many volume snippets we have
    for (size_t vsnc=0; vsnc<pdd->number_to_volumesnippet.size(); vsnc++) // loop over volume snippet numbers on this level
      {
	std::set<size_t> globaldofs_in_volumesnippet; // empty set
	for (auto vsnf : pdd->volumesnippet_coarse_to_fine[vsnc]) // loop over all fine volume snippet numbers contained in the coarse volume snippet
	  for (auto p : pdd->fine_pdd->number_to_volumesnippet[vsnf]) // loop over all subdomains on fine level the volume snippet is in
	    for (size_t k=0; k<fine_pdddm->n_basis_vectors_per_subdomain[p]; k++) // loop over all basis vectors in this subdomain
	      globaldofs_in_volumesnippet.insert(fine_pdddm->basisvector_local_to_global[p][k]); // now we have fond a global dof
	// now we know all global dofs for this volume snippet
	for (auto g : globaldofs_in_volumesnippet)
	  volumesnippetnumber_local_to_global[vsnc].push_back(g); // enter all global dof in local to global map
      }
    // now that we have local to global build up global to local
    volumesnippetnumber_global_to_local.resize(n_global_dofs);
    for (size_t vsn=0; vsn<volumesnippetnumber_local_to_global.size(); vsn++) // loop over volume snippet numbers
      for (size_t k=0; k<volumesnippetnumber_local_to_global[vsn].size(); k++) // loop over all dofs in that volume snippet
	volumesnippetnumber_global_to_local[volumesnippetnumber_local_to_global[vsn][k]][vsn] = k;

    // There are also skeleton snippets on the fine level which are inside volume snippets on this level but they will not add any dofs

    // now for the skeleton snippets
    skeletonsnippetnumber_local_to_global.resize(pdd->number_to_skeletonsnippet.size()); // resize to correct length
    for (size_t ssnc=0; ssnc<pdd->number_to_skeletonsnippet.size(); ssnc++) // loop over skeleton snippet numbers on this level
      {
	std::set<size_t> globaldofs_in_skeletonsnippet; // empty set
	for (auto ssnf : pdd->skeletonsnippet_coarse_to_fine[ssnc]) // loop over all fine skeleton snippet numbers contained in this coarse skeleton snippet
	  {
	    auto& skeletonsnippet = pdd->fine_pdd->number_to_skeletonsnippet[ssnf]; // fine skeleton snippet with number ssnf
	    //std::cout << pdd->vs_to_str(skeletonsnippet) << std::endl;
	    auto it = skeletonsnippet.begin(); // points to first v snippet number
	    size_t left_volumesnippetnumber = *it;
	    ++it;
	    size_t right_volumesnippetnumber = *it; // so now we have the numbers of the two volume snippets to left and right
	    auto& left_volumesnippet = pdd->fine_pdd->number_to_volumesnippet[left_volumesnippetnumber];
	    auto& right_volumesnippet = pdd->fine_pdd->number_to_volumesnippet[right_volumesnippetnumber];
	    for (auto p : left_volumesnippet) // loop over all subdomains on fine level the left volume snippet is in
	      for (size_t k=0; k<fine_pdddm->n_basis_vectors_per_subdomain[p]; k++) // loop over all basis vectors in this fine subdomain
		globaldofs_in_skeletonsnippet.insert(fine_pdddm->basisvector_local_to_global[p][k]); // now we have fond a global dof
	    for (auto p : right_volumesnippet) // loop over all subdomains on fine level the right volume snippet is in
	      for (size_t k=0; k<fine_pdddm->n_basis_vectors_per_subdomain[p]; k++) // loop over all basis vectors in this subdomain
		globaldofs_in_skeletonsnippet.insert(fine_pdddm->basisvector_local_to_global[p][k]); // now we have fond a global dof
	  }
	// now we know all global dofs for this skeleton snippet
	for (auto g : globaldofs_in_skeletonsnippet)
	  skeletonsnippetnumber_local_to_global[ssnc].push_back(g); // enter all global dof in local to global map
      }
    // now that we have local to global build up global to local
    skeletonsnippetnumber_global_to_local.resize(n_global_dofs);
    for (size_t ssn=0; ssn<skeletonsnippetnumber_local_to_global.size(); ssn++) // loop over volume snippet numbers
      for (size_t k=0; k<skeletonsnippetnumber_local_to_global[ssn].size(); k++) // loop over all dofs in that volume snippet
	skeletonsnippetnumber_global_to_local[skeletonsnippetnumber_local_to_global[ssn][k]][ssn] = k;

    // now subdomain local <-> global for this level
    size_t n_subdomains = pdd->subdomains;
    subdomainnumber_local_to_global.resize(n_subdomains); // resize for fill-in
    subdomainnumber_global_to_local.resize(n_global_dofs); // n_global_dofs is still the same ...
    for (size_t p=0; p<n_subdomains; ++p) // loop over all subdomains
      {
	std::set<size_t> globaldofs; // find global dofs for this subdomain
	
	for (auto volumesnippetnumber : pdd->subdomain_to_volumesnippetnumber[p]) // loop over all volumesnippets in that subdomain
	  for (auto g : volumesnippetnumber_local_to_global[volumesnippetnumber]) // loop over all global dof numbers in that snippet
	    globaldofs.insert(g); // add global dof to subdomain

	// now we can build up local to global map
	for (auto g : globaldofs)
	  subdomainnumber_local_to_global[p].push_back(g);

	// and global to local
	for (size_t i=0; i<subdomainnumber_local_to_global[p].size(); i++)
	  subdomainnumber_global_to_local[subdomainnumber_local_to_global[p][i]][p] = i;
      }
    // boundary information
    subdomain_boundary.resize(n_subdomains);
    for (size_t p=0; p<n_subdomains; ++p) // loop over all coarse subdomains
      {
	subdomain_boundary[p].resize(subdomainnumber_local_to_global[p].size(),false);
	for (size_t i=0; i<subdomainnumber_local_to_global[p].size(); ++i) // loop over all local dofs in this subdomain
	  {
	    auto g = subdomainnumber_local_to_global[p][i]; // the global dof
	    auto q = fine_pdddm->basisvector_global_to_subdomainnumber[g]; // the fine level subdomain this dof originates from
	    if (pdd->subdomain_coarse_to_fine_map[p].count(q)==0) // is this fine subdomain mapped to the coarse subdomain?
	      subdomain_boundary[p][i] = true;
	  }
      }

    // construct subdomain graph on this level
    p_subdomain_graph = std::shared_ptr<WeightedGraph>(new WeightedGraph()); // allocate graph
    WeightedGraph& g=*p_subdomain_graph; // for simplicity
    g.graph.resize(pdd->subdomains);
    g.vertex_weights.resize(pdd->subdomains,0);
    g.edge_weights.resize(pdd->subdomains);
    for (size_t j=0; j<subdomainnumber_global_to_local.size(); j++) // loop over global dof numbers
      for (auto itp = subdomainnumber_global_to_local[j].begin(); itp!=subdomainnumber_global_to_local[j].end(); ++itp) // loop over subdomain numbers sharing this degree of freedom
	{
	  g.vertex_weights[itp->first]++; // increment vertex weight
	  for (auto itq = subdomainnumber_global_to_local[j].begin(); itq!=subdomainnumber_global_to_local[j].end(); ++itq) // loop over subdomain numbers sharing this degree of freedom
	    {
	      g.graph[itp->first].insert(itq->first); // couple the two subdomains because they share a degree of freedom
	      if (itp->first!=itq->first) // edge weights are only calculated for offdiagonal edges
		if (g.edge_weights[itp->first].find(itq->first)==g.edge_weights[itp->first].end())
		  g.edge_weights[itp->first][itq->first] = 1; // new entry
		else
		  g.edge_weights[itp->first][itq->first] += 1; // increment existing entry
	    }
	}

    // check the subdomain graph, if we have the same connectivity
    if (pdd->subdomains!=p_subdomain_graph->graph.size())
      DUNE_THROW(Dune::Exception, "DomainDecompositionDOFMapper: size mismatch -X");
    std::vector<std::set<size_t>> mygraph(pdd->subdomains);
    for (size_t p=0; p<pdd->subdomains; ++p)
      for (auto g : subdomainnumber_local_to_global[p])
	for (auto& entry : subdomainnumber_global_to_local[g])
	  mygraph[p].insert(entry.first);
    for (size_t p=0; p<pdd->subdomains; ++p)
      if (mygraph[p]!=p_subdomain_graph->graph[p])
	{
	  std::cout << p << "  mygraph"; for (auto q : mygraph[p]) std::cout << " " << q; std::cout << std::endl; 
	  std::cout << p << " pddgraph"; for (auto q : p_subdomain_graph->graph[p]) std::cout << " " << q; std::cout << std::endl; 
	  DUNE_THROW(Dune::Exception, "DomainDecompositionDOFMapper: subdomain graphs are not the same");
	}
  }
  
  void print_info (int verbosity=0)
  {
    if (verbosity==0) return;
    std::cout << "DomainDecompositionDOFMapper [" << std::endl;
    auto& dd = *pdd; // use reference, its nicer
    if (verbosity>=1)
      std::cout << n_global_dofs << " global dofs with blocksize " << blocksize << std::endl;
    if (verbosity>=2)
      for (size_t i=0; i<volumesnippetnumber_local_to_global.size(); i++)
	{
	  std::cout << "volume snippet number " << i << "=(";
	  for (auto p : dd.number_to_volumesnippet[i]) std::cout << " " << p;
	  std::cout << " ) has " << volumesnippetnumber_local_to_global[i].size() << " degrees of freedom (blocks)";
	  if (verbosity>=2) for (auto g : volumesnippetnumber_local_to_global[i]) std::cout << " " << g;
	  std::cout << std::endl;
	}
    if (verbosity>=3)
      for (size_t g=0; g<volumesnippetnumber_global_to_local.size(); g++)
        {
          std::cout << "global index " << g << ": ";
          for (auto entry : volumesnippetnumber_global_to_local[g])
            std::cout << " (" << entry.first << "|" << entry.second <<")";
          std::cout << std::endl;
        }
    if (verbosity>=2)
      for (size_t i=0; i<skeletonsnippetnumber_local_to_global.size(); i++)
	{
	  std::cout << "skeleton snippet number " << i << "=(";
	  for (auto p : dd.number_to_skeletonsnippet[i]) std::cout << " " << p;
	  std::cout << " ) has " << skeletonsnippetnumber_local_to_global[i].size() << " degrees of freedom (blocks)";
	  if (verbosity>=2) for (auto g : skeletonsnippetnumber_local_to_global[i]) std::cout << " " << g;
	  std::cout << std::endl;
	}
    if (verbosity>=2)
      for (size_t vsn=0; vsn<volumesnippetnumber_local_to_global.size(); vsn++) // loop over volume snippet numbers
	std::cout << "volume snippet " << vsn << " has " << volumesnippetnumber_local_to_global[vsn].size() << " dofs" << std::endl;
    if (verbosity>=2)
      for (size_t ssn=0; ssn<skeletonsnippetnumber_local_to_global.size(); ssn++) // loop over volume snippet numbers
	std::cout << "skeleton snippet " << ssn << " has " << skeletonsnippetnumber_local_to_global[ssn].size() << " dofs" << std::endl;
    if (verbosity>=2)
      for (size_t p=0; p<subdomainnumber_local_to_global.size(); ++p)
	{
	  std::cout << "subdomain " << p << " has " << subdomainnumber_local_to_global[p].size() << " dofs" << std::endl;
	  size_t n=0;
	  if (subdomain_boundary.size()==subdomainnumber_local_to_global.size())
	    {
	      std::cout << "   interior dofs are:";
	      for (size_t i=0; i<subdomainnumber_local_to_global[p].size(); ++i) // loop over all local dofs in this subdomain
		if (!subdomain_boundary[p][i])
		  {
		    std::cout << " " << subdomainnumber_local_to_global[p][i];
		    n++;
		  }
	      std::cout << " (" << n << ")" << std::endl;
	      n=0;
	      std::cout << "   boundary dofs are:";
	      for (size_t i=0; i<subdomainnumber_local_to_global[p].size(); ++i) // loop over all local dofs in this subdomain
		if (subdomain_boundary[p][i])
		  {
		    std::cout << " " << subdomainnumber_local_to_global[p][i];
		    n++;
		  }
	      std::cout << " (" << n << ")" << std::endl;
	    }
	}
    std::cout << "]" << std::endl;
  }

};


// set up sparse matrices for subsets of elements given by a local<->global mapping
// using the global stiffness matrix as a template
template<typename ISTLM>
std::vector<std::shared_ptr<ISTLM>> set_up_local_matrices (const ISTLM& Aglobal,
							   const std::vector<std::vector<size_t>>& local_to_global,
							   const std::vector<std::map<size_t,size_t>>& global_to_local)
{
  auto n_pieces = local_to_global.size();
  int avg = Aglobal.nonzeroes()/Aglobal.N()+2;
  std::vector<std::shared_ptr<ISTLM>> matrices(n_pieces); // for each volumesnippet we wil have one matrix

  for (size_t p=0; p<n_pieces; ++p) // loop over all pieces, for each we will have a matrix 
    {
      auto n = local_to_global[p].size(); // number of DOFs in piece p
      //std::cout << "building matrix in piece " << p << " size " << n << " times " << n << " avg row size " << avg << std::endl;
      auto pA = std::shared_ptr<ISTLM>(new ISTLM(n,n,avg,0.15,ISTLM::implicit)); // a new sparse matrix
      for (size_t i=0; i<n; i++)
        {
          auto iglobal = local_to_global[p][i]; // row in global matrix
          auto cIt = Aglobal[iglobal].begin();
          auto cEndIt = Aglobal[iglobal].end();
          for (; cIt!=cEndIt; ++cIt) // loop over row in global matrix
            {
              auto jglobal = cIt.index();
              auto it = global_to_local[jglobal].find(p); // ask for global column index being in snippet s
              if (it!=global_to_local[jglobal].end())
                {
                  // ok, now we have an entry
                  auto j = it->second; // global to local index
                  pA->entry(i,j) = 0.0;
                }
            }
        }
      auto stats = pA->compress();
      matrices[p] = pA; // store pointer
    }

  // return all the matrices
  return matrices;
}

// set up sparse matrices for subsets of elements given by a local<->global mapping
// using the global stiffness matrix as a template
template<typename ISTLM>
std::shared_ptr<ISTLM> set_up_local_matrix (const ISTLM& Aglobal, // global stiffness matrix
					    size_t p, // the piece to assemble
					    const std::vector<std::vector<size_t>>& local_to_global,
					    const std::vector<std::map<size_t,size_t>>& global_to_local)
{
  int avg = Aglobal.nonzeroes()/Aglobal.N()+2;

  auto n = local_to_global[p].size(); // number of DOFs in piece p
  //std::cout << "building matrix in piece " << p << " size " << n << " times " << n << " avg row size " << avg << std::endl;
  auto pA = std::shared_ptr<ISTLM>(new ISTLM(n,n,avg,0.15,ISTLM::implicit)); // a new sparse matrix
  for (size_t i=0; i<n; i++)
    {
      auto iglobal = local_to_global[p][i]; // row in global matrix
      auto cIt = Aglobal[iglobal].begin();
      auto cEndIt = Aglobal[iglobal].end();
      for (; cIt!=cEndIt; ++cIt) // loop over row in global matrix
	{
	  auto jglobal = cIt.index();
	  auto it = global_to_local[jglobal].find(p); // ask for global column index being in snippet s
	  if (it!=global_to_local[jglobal].end())
	    {
	      // ok, now we have an entry
	      auto j = it->second; // global to local index
	      pA->entry(i,j) = 0.0;
	    }
	}
    }
  auto stats = pA->compress();

  // return the matrix
  return pA;
}


// set up sparse matrices for subsets of elements given by local<->global mapping
// using the global stiffness matrix as a template
template<typename ISTLM, typename C1>
void assemble_snippets_to_subdomain_matrices (const std::vector<C1>& subdomain_to_snippetnumber,
					      const std::vector<std::vector<size_t>>& local_to_global,
					      const std::vector<std::map<size_t,size_t>>& global_to_local,
					      std::vector<std::shared_ptr<ISTLM>>& snippet_matrices,
					      std::vector<std::shared_ptr<ISTLM>>& subdomain_matrices) // this is the result to be filled
{
  // we assume that the subdomain matrices are set up and zeroed
  for (size_t p=0; p<subdomain_matrices.size(); p++) // loop over subdomains
    {
      auto& subdomainA = *subdomain_matrices[p]; // get matrix in subdomain p
      for (auto snippetnumber : subdomain_to_snippetnumber[p]) // loop over all snippets making up subdomain p
	{
	  auto& snippetA = *snippet_matrices[snippetnumber]; // get the matrix for this snippet
	  // now we want to scatter that snippet matrix to the subdomain matrix
	  for (auto rIt=snippetA.begin(); rIt!=snippetA.end(); ++rIt)
	    {
	      auto isnippet = rIt.index();
	      auto iglobal = local_to_global[snippetnumber][isnippet];
	      auto iti = global_to_local[iglobal].find(p);
	      auto isubdomain = iti->second;
	      auto cIt = rIt->begin();
	      auto cEndIt = rIt->end();
	      for (; cIt!=cEndIt; ++cIt)
		{
		  auto jsnippet = cIt.index();
		  auto jglobal = local_to_global[snippetnumber][jsnippet];
		  auto itj = global_to_local[jglobal].find(p);
		  auto jsubdomain = itj->second;
		  subdomainA[isubdomain][jsubdomain] += (*cIt);
		}
	    }
	}
    }
}


// set up sparse matrices for subsets of elements given by local<->global mapping
// using the global stiffness matrix as a template
template<typename ISTLM, typename C>
void assemble_snippets_to_subdomain_matrix (size_t p, // the subdomain to assemble to
					    const std::vector<C>& subdomain_to_snippetnumber,
					    const std::vector<std::vector<size_t>>& local_to_global,
					    const std::vector<std::map<size_t,size_t>>& global_to_local,
					    std::vector<std::shared_ptr<ISTLM>>& snippet_matrices,
					    std::shared_ptr<ISTLM> subdomain_matrix
					    ) // this is the result to be filled
{
  // we assume that the subdomain matrix exists and is appropriately initialized
  auto& subdomainA = *subdomain_matrix; // get matrix in subdomain p
  for (auto snippetnumber : subdomain_to_snippetnumber[p]) // loop over all snippets making up subdomain p
    {
      auto& snippetA = *snippet_matrices[snippetnumber]; // get the matrix for this snippet
      // now we want to scatter that snippet matrix to the subdomain matrix
      for (auto rIt=snippetA.begin(); rIt!=snippetA.end(); ++rIt)
	{
	  auto isnippet = rIt.index();
	  auto iglobal = local_to_global[snippetnumber][isnippet];
	  auto iti = global_to_local[iglobal].find(p);
	  if (iti==global_to_local[iglobal].end())
	    {
	      std::cout << "global index " << iglobal << " not found in subdomain " << p << std::endl;
	      DUNE_THROW(Dune::Exception, "access error.");
	    }
	  auto isubdomain = iti->second;
	  auto cIt = rIt->begin();
	  auto cEndIt = rIt->end();
	  for (; cIt!=cEndIt; ++cIt)
	    {
	      auto jsnippet = cIt.index();
	      auto jglobal = local_to_global[snippetnumber][jsnippet];
	      auto itj = global_to_local[jglobal].find(p);
	      if (iti==global_to_local[iglobal].end())
		{
		  std::cout << "global index " << jglobal << " not found in subdomain " << p << std::endl;
		  DUNE_THROW(Dune::Exception, "access error.");
		}
	      auto jsubdomain = itj->second;
	      if (!subdomainA.exists(isubdomain,jsubdomain))
		{
		  std::cout << "entry " <<  isubdomain << "," << jsubdomain << " not in matrix in subdomain " << p
			    << " global: " << iglobal << "," << jglobal << " value=" << *cIt
			    << " snippetnumber=" << snippetnumber
			    << std::endl;
		  if (std::abs((*cIt)[0][0])>1e-12)
		    DUNE_THROW(Dune::Exception, "access error.");
		}
	      else
		subdomainA[isubdomain][jsubdomain] += (*cIt);
	    }
	}
    }
}


// a variant that omits the interior elements; for volume snippets
template<typename ISTLM, typename C>
void assemble_snippets_to_subdomain_matrix_geneo_volume (size_t p, // the subdomain to assemble to
					    const std::vector<C>& subdomain_to_snippetnumber,
					    const std::vector<std::set<size_t>>& number_to_volumesnippet,
					    const std::vector<std::vector<size_t>>& local_to_global,
					    const std::vector<std::map<size_t,size_t>>& global_to_local,
					    std::vector<std::shared_ptr<ISTLM>>& snippet_matrices,
					    ISTLM& subdomain_matrix
					    ) // this is the result to be filled
{
  // we assume that the subdomain matrix exists and is appropriately initialized
  auto& subdomainA = subdomain_matrix; // get matrix in subdomain p
  for (auto snippetnumber : subdomain_to_snippetnumber[p]) // loop over all snippets making up subdomain p
    {
      auto& snippetA = *snippet_matrices[snippetnumber]; // get the matrix for this snippet
      // now we want to scatter that snippet matrix to the subdomain matrix
      auto snippet = number_to_volumesnippet[snippetnumber];
      if (snippet.size()==1 && snippet.count(p)==1)
	continue; // omit assembling the interior!
      for (auto rIt=snippetA.begin(); rIt!=snippetA.end(); ++rIt)
	{
	  auto isnippet = rIt.index();
	  auto iglobal = local_to_global[snippetnumber][isnippet];
	  auto iti = global_to_local[iglobal].find(p);
	  if (iti==global_to_local[iglobal].end())
	    {
	      std::cout << "global index " << iglobal << " not found in subdomain " << p << std::endl;
	      DUNE_THROW(Dune::Exception, "access error.");
	    }
	  auto isubdomain = iti->second;
	  auto cIt = rIt->begin();
	  auto cEndIt = rIt->end();
	  for (; cIt!=cEndIt; ++cIt)
	    {
	      auto jsnippet = cIt.index();
	      auto jglobal = local_to_global[snippetnumber][jsnippet];
	      auto itj = global_to_local[jglobal].find(p);
	      if (iti==global_to_local[iglobal].end())
		{
		  std::cout << "global index " << jglobal << " not found in subdomain " << p << std::endl;
		  DUNE_THROW(Dune::Exception, "access error.");
		}
	      auto jsubdomain = itj->second;
	      if (!subdomainA.exists(isubdomain,jsubdomain))
		{
		  std::cout << "entry " <<  isubdomain << "," << jsubdomain << " not in matrix in subdomain " << p
			    << " global: " << iglobal << "," << jglobal << " value=" << *cIt
			    << " snippetnumber=" << snippetnumber
			    << std::endl;
		  if (std::abs((*cIt)[0][0])>1e-12)
		    DUNE_THROW(Dune::Exception, "access error.");
		}
	      else
		subdomainA[isubdomain][jsubdomain] += (*cIt);
	    }
	}
    }
}

// a variant that omits the faces touching the interior; for skeleton snippets
template<typename ISTLM, typename C>
void assemble_snippets_to_subdomain_matrix_geneo_skeleton (size_t p, // the subdomain to assemble to
					    const std::vector<C>& subdomain_to_snippetnumber,
					    const std::vector<std::set<size_t>>& number_to_volumesnippet,
					    const std::vector<std::set<size_t>>& number_to_skeletonsnippet,
					    const std::vector<std::vector<size_t>>& local_to_global,
					    const std::vector<std::map<size_t,size_t>>& global_to_local,
					    std::vector<std::shared_ptr<ISTLM>>& snippet_matrices,
					    ISTLM& subdomain_matrix
					    ) // this is the result to be filled
{
  // we assume that the subdomain matrix exists and is appropriately initialized
  auto& subdomainA = subdomain_matrix; // get matrix in subdomain p
  for (auto snippetnumber : subdomain_to_snippetnumber[p]) // loop over all snippets making up subdomain p
    {
      auto& snippetA = *snippet_matrices[snippetnumber]; // get the matrix for this snippet
      // now we want to scatter that snippet matrix to the subdomain matrix

      auto skeletonsnippet = number_to_skeletonsnippet[snippetnumber];
      auto it = skeletonsnippet.begin(); // points to first snippet number
      size_t left_volumesnippetnumber = *it;
      ++it;
      size_t right_volumesnippetnumber = *it; // so now we have the numbers of the two volume snippets to left and right
      auto& left_volumesnippet = number_to_volumesnippet[left_volumesnippetnumber];
      auto& right_volumesnippet = number_to_volumesnippet[right_volumesnippetnumber];

      // skip, if the interior is on either side of the skeleton snippet
      if (left_volumesnippet.size()==1 && left_volumesnippet.count(p)==1)
	continue; // omit assembling the interior!
      if (right_volumesnippet.size()==1 && right_volumesnippet.count(p)==1)
	continue; // omit assembling the interior!

      for (auto rIt=snippetA.begin(); rIt!=snippetA.end(); ++rIt)
	{
	  auto isnippet = rIt.index();
	  auto iglobal = local_to_global[snippetnumber][isnippet];
	  auto iti = global_to_local[iglobal].find(p);
	  if (iti==global_to_local[iglobal].end())
	    {
	      std::cout << "global index " << iglobal << " not found in subdomain " << p << std::endl;
	      DUNE_THROW(Dune::Exception, "access error.");
	    }
	  auto isubdomain = iti->second;
	  auto cIt = rIt->begin();
	  auto cEndIt = rIt->end();
	  for (; cIt!=cEndIt; ++cIt)
	    {
	      auto jsnippet = cIt.index();
	      auto jglobal = local_to_global[snippetnumber][jsnippet];
	      auto itj = global_to_local[jglobal].find(p);
	      if (iti==global_to_local[iglobal].end())
		{
		  std::cout << "global index " << jglobal << " not found in subdomain " << p << std::endl;
		  DUNE_THROW(Dune::Exception, "access error.");
		}
	      auto jsubdomain = itj->second;
	      if (!subdomainA.exists(isubdomain,jsubdomain))
		{
		  std::cout << "entry " <<  isubdomain << "," << jsubdomain << " not in matrix in subdomain " << p
			    << " global: " << iglobal << "," << jglobal << " value=" << *cIt
			    << " snippetnumber=" << snippetnumber
			    << std::endl;
		  if (std::abs((*cIt)[0][0])>1e-12)
		    DUNE_THROW(Dune::Exception, "access error.");
		}
	      else
		subdomainA[isubdomain][jsubdomain] += (*cIt);
	    }
	}
    }
}

#endif // Udune_ftworth_HH
