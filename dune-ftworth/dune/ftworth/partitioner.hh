#ifndef Udune_ftworth_partitioner_HH
#define Udune_ftworth_partitioner_HH

#if HAVE_PARMETIS

#include <parmetis.h>

#include <algorithm>
#include <vector>
#include <dune/grid/utility/globalindexset.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>

// only enable for ParMETIS because the implementation uses functions that
// are not emulated by scotch
#ifdef PARMETIS_MAJOR_VERSION

/** \brief Returns a vector with a partition number for each element
 *
 * \author Benjamin Bykowski and adapted by Peter Bastian
 *
 * \param gv The grid view to be partitioned
 * \param mpihelper The MPIHelper object, needed to get the MPI communicator. This is needed by the function ParMETIS_V3_PartMeshKway and can unfortunately not be omitted
 * \param parts number of subdomains desired
 *
 * \return std::vector with one uint per All_Partition element.  For each element, the entry is the
 *    number of the partition the element is assigned to. This number is greater or equal zero and smaller as parts. 
 *    No partitioning is done, only this vector is computed
 */
template<class GridView>
std::vector<size_t> parmetis_partitioning (const GridView& gv, const Dune::MPIHelper& mpihelper, int parts) {

#if PARMETIS_MAJOR_VERSION > 3
  typedef idx_t idx_type;
  typedef ::real_t real_type;
#else
  typedef int idx_type;
  typedef float real_type;
#endif // PARMETIS_MAJOR_VERSION > 3

  const unsigned numElements = gv.size(0);

  std::vector<idx_type> part(numElements);

  // Setup parameters for ParMETIS
  idx_type wgtflag = 0;                                  // we don't use weights
  idx_type numflag = 0;                                  // we are using C-style arrays
  idx_type ncon = 1;                                     // number of balance constraints
  idx_type ncommonnodes = 2;                             // number of nodes elements must have in common to be considered adjacent to each other
  idx_type options[4] = {0, 0, 0, 0};                    // use default values for random seed, output and coupling
  idx_type edgecut;                                      // will store number of edges cut by partition
  idx_type nparts = parts;                               // number of partitions to create is a parameter
  std::vector<real_type> tpwgts(ncon*nparts, 1./nparts); // load per subdomain and weight (same load on every process)
  std::vector<real_type> ubvec(ncon, 1.05);              // weight tolerance (same weight tolerance for every weight there is)

  // The difference elmdist[i+1] - elmdist[i] is the number of nodes that are on process i
  std::vector<idx_type> elmdist(nparts+1);
  elmdist[0] = 0;
  std::fill(elmdist.begin()+1, elmdist.end(), gv.size(0)); // all elements are on process zero

  // Create and fill arrays "eptr", where eptr[i] is the number of vertices that belong to the i-th element, and
  // "eind" contains the vertex-numbers of the i-the element in eind[eptr[i]] to eind[eptr[i+1]-1]
  std::vector<idx_type> eptr, eind;
  int numVertices = 0;
  eptr.push_back(numVertices);

  for (const auto& element : elements(gv, Dune::Partitions::interior)) {
    const size_t curNumVertices = Dune::referenceElement<double,GridView::dimension>(element.type()).size(GridView::dimension);

    numVertices += curNumVertices;
    eptr.push_back(numVertices);

    for (size_t k = 0; k < curNumVertices; ++k)
      eind.push_back(gv.indexSet().subIndex(element, k, GridView::dimension));
  }

  // Partition mesh using ParMETIS
  if (0 == mpihelper.rank()) {
    MPI_Comm comm = Dune::MPIHelper::getLocalCommunicator();

#if PARMETIS_MAJOR_VERSION >= 4
    const int OK =
#endif
      ParMETIS_V3_PartMeshKway(elmdist.data(), eptr.data(), eind.data(), NULL, &wgtflag, &numflag,
                               &ncon, &ncommonnodes, &nparts, tpwgts.data(), ubvec.data(),
                               options, &edgecut, part.data(), &comm);

#if PARMETIS_MAJOR_VERSION >= 4
    if (OK != METIS_OK)
      DUNE_THROW(Dune::Exception, "ParMETIS returned an error code.");
#endif
  }
  std::vector<size_t> retpart(part.size());
  for (size_t i=0; i<part.size(); ++i) retpart[i] = part[i];
  return retpart;
}



/** \brief Returns a vector with a partition number for each element
 *
 * \author Benjamin Bykowski and adapted by Peter Bastian
 *
 * \param gv The grid view to be partitioned
 * \param comm an MPI communicator. This is needed by the function ParMETIS_V3_PartMeshKway and can unfortunately not be omitted
 * \param parts number of subdomains desired
 *
 * \return std::vector with one uint per All_Partition element.  For each element, the entry is the
 *    number of the partition the element is assigned to. This number is greater or equal zero and smaller as parts. 
 *    No partitioning is done, only this vector is computed
 */
template<class GridView>
std::vector<size_t> parmetis_partitioning (const GridView& gv, MPI_Comm comm, int parts) {

#if PARMETIS_MAJOR_VERSION > 3
  typedef idx_t idx_type;
  typedef ::real_t real_type;
#else
  typedef int idx_type;
  typedef float real_type;
#endif // PARMETIS_MAJOR_VERSION > 3

  const unsigned numElements = gv.size(0);

  std::vector<idx_type> part(numElements);

  // Setup parameters for ParMETIS
  idx_type wgtflag = 0;                                  // we don't use weights
  idx_type numflag = 0;                                  // we are using C-style arrays
  idx_type ncon = 1;                                     // number of balance constraints
  idx_type ncommonnodes = 2;                             // number of nodes elements must have in common to be considered adjacent to each other
  idx_type options[4] = {0, 0, 0, 0};                    // use default values for random seed, output and coupling
  idx_type edgecut;                                      // will store number of edges cut by partition
  idx_type nparts = parts;                               // number of partitions to create is a parameter
  std::vector<real_type> tpwgts(ncon*nparts, 1./nparts); // load per subdomain and weight (same load on every process)
  std::vector<real_type> ubvec(ncon, 1.05);              // weight tolerance (same weight tolerance for every weight there is)

  // The difference elmdist[i+1] - elmdist[i] is the number of nodes that are on process i
  std::vector<idx_type> elmdist(nparts+1);
  elmdist[0] = 0;
  std::fill(elmdist.begin()+1, elmdist.end(), gv.size(0)); // all elements are on process zero

  // Create and fill arrays "eptr", where eptr[i] is the number of vertices that belong to the i-th element, and
  // "eind" contains the vertex-numbers of the i-the element in eind[eptr[i]] to eind[eptr[i+1]-1]
  std::vector<idx_type> eptr, eind;
  int numVertices = 0;
  eptr.push_back(numVertices);

  for (const auto& element : elements(gv, Dune::Partitions::interior)) {
    const size_t curNumVertices = Dune::referenceElement<double,GridView::dimension>(element.type()).size(GridView::dimension);

    numVertices += curNumVertices;
    eptr.push_back(numVertices);

    for (size_t k = 0; k < curNumVertices; ++k)
      eind.push_back(gv.indexSet().subIndex(element, k, GridView::dimension));
  }

  // Partition mesh using ParMETIS
  int rank;
  MPI_Comm_rank (comm, &rank);
  if (rank==0) {

#if PARMETIS_MAJOR_VERSION >= 4
    const int OK =
#endif
      ParMETIS_V3_PartMeshKway(elmdist.data(), eptr.data(), eind.data(), NULL, &wgtflag, &numflag,
                               &ncon, &ncommonnodes, &nparts, tpwgts.data(), ubvec.data(),
                               options, &edgecut, part.data(), &comm);

#if PARMETIS_MAJOR_VERSION >= 4
    if (OK != METIS_OK)
      DUNE_THROW(Dune::Exception, "ParMETIS returned an error code.");
#endif
  }
  std::vector<size_t> retpart(part.size());
  for (size_t i=0; i<part.size(); ++i) retpart[i] = part[i];
  return retpart;
}


/** \brief Returns a vector with a partition number for each element
 *
 * Compute bounding box and subdivide the bounding box into the given number of parts.
 * Then assign elements to these parts according to their center.
 *
 * \author Peter Bastian
 *
 * \param gv The grid view to be partitioned
 * \param parts number of subdomains in each direction
 *
 * \return std::vector with one uint per All_Partition element.  For each element, the entry is the
 *    number of the partition the element is assigned to. This number is greater or equal zero and smaller as parts. 
 *    No partitioning is done, only this vector is computed
 */
template<class GridView>
std::vector<size_t> coordinate_partitioning (const GridView& gv, std::vector<size_t> parts)
{
    // grid information and procs
    const size_t dim = GridView::dimension;
    auto& indexset = gv.indexSet(); // to attach data to elements
    if (parts.size()!=dim)
      DUNE_THROW(Dune::Exception, "coordinate_partitioning: parts.size()!=dim");
    size_t procs = 1; for (auto i : parts) procs *= i;

    // determine bounding box of the domain
    using ctype = typename GridView::Grid::ctype;
    std::vector<ctype> cmin(dim,1e100);
    std::vector<ctype> cmax(dim,-1e100);
    for (const auto& vertex : vertices(gv))
      for (int i=0; i<dim; i++)
	{
	  cmin[i] = std::min(cmin[i],vertex.geometry().corner(0)[i]);
	  cmax[i] = std::max(cmax[i],vertex.geometry().corner(0)[i]);
	}
    std::vector<ctype> H(dim);
    for (int i=0; i<dim; i++) H[i] = (cmax[i]-cmin[i])/parts[i];

    // assign elements
    std::vector<size_t> partitioning(gv.size(0));
    std::vector<size_t> N(dim,0);
    N[0]=1; for (int i=1; i<dim; i++) N[i] = N[i-1]*parts[i-1];
    for (const auto& element : elements(gv))
      {
	auto center = element.geometry().center();
	std::vector<size_t> pos(dim);
	for (int i=0; i<dim; i++) pos[i] = floor((center[i]-cmin[i])/H[i]);
	for (int i=0; i<dim; i++)
	  if (pos[i]<0 || pos[i]>=parts[i])
	    {
	      std::cout << "i=" << i << " pos=" << pos[i] << std::endl;
	      DUNE_THROW(Dune::Exception, "coordinate_partitioning: pos out of range");
	    }    
	size_t part = 0;
	for (int i=0; i<dim; i++) part += pos[i]*N[i];
	if (part<0 || part>=procs)
	  {
	    std::cout << "center=" << center << std::endl;
	    std::cout << "H="; for (auto h : H) std::cout << h << " "; std::cout << std::endl;
	    std::cout << "pos="; for (auto i : pos) std::cout << i << " "; std::cout << std::endl;
	    std::cout << "part=" << part << std::endl;
	    DUNE_THROW(Dune::Exception, "coordinate_partitioning: part out of range");
	  }
	partitioning[indexset.index(element)] = part;
      }

  return partitioning;
}


/**
 * \brief partition a graph representing subdomains
 *
 * \param subdomain_graph graph represented by std::vector<std::set<size_t>>; a vertex may be coupled to itself
 * \param vertex_weights weight for each vertex in the input graph
 * \param edge_weights weight for each edge in the input graph, only non diagonal entries are used
 * \param mpihelper is required by parmetis
 * \param parts number of partitions in the output
 *
 * Partition the input graph using the method ParMETIS_V3_PartKway of ParMetis. Weights must be given because
 * otherwise the partitioning might be lousy.
 *
 * \return std::vector<unsigned>> for each element stores a set of subomains containing this element
 */
std::vector<size_t> parmetis_sequential_graph_partitioning (const std::vector<std::set<size_t>>& subdomain_graph,
							      const std::vector<size_t> vertex_weights, // for each entry in subdomain_graph
							      const std::vector<std::map<size_t,size_t>> edge_weights, // for each non-diagonal entry in subdomain_graph
							      const Dune::MPIHelper& mpihelper,
							      int parts)
{
#if PARMETIS_MAJOR_VERSION > 3
  typedef idx_t idx_type;
  typedef ::real_t real_type;
#else
  typedef int idx_type;
  typedef float real_type;
#endif // PARMETIS_MAJOR_VERSION > 3

  // compute some quantities
  size_t nprocs = 1; // number of processors the graph is already distributed on
  size_t n = subdomain_graph.size(); // the number of vertices in the input graph;
  size_t m = 0; // the total number of edges in the input graph
  for (unsigned p=0; p<n; p++)
    for (auto q : subdomain_graph[p])
      if (q!=p)
	m++;

  // prepare ParMetis input
  std::vector<idx_type> vtxdist(nprocs+1); // vertex distribution, everything is sequential here
  vtxdist[0] = 0;
  vtxdist[1] = n;
  std::vector<idx_type> xadj(n+1); // CRS format
  std::vector<idx_type> vwgt(n);
  std::vector<idx_type> adjncy(m);
  std::vector<idx_type> adjwgt(m);
  size_t count=0;
  for (unsigned p=0; p<n; p++)
    {
      xadj[p] = count;
      vwgt[p] = vertex_weights[p];
      for (auto q : subdomain_graph[p])
	if (q!=p)
	  {
	    adjncy[count] = q;
	    adjwgt[count] = edge_weights[p].find(q)->second;
	    count++;
	  }
    }
  xadj[n] = count;
  idx_type wgtflag = 0; // we use weights for vertices and edges
  idx_type numflag = 0; // we are using C-style arrays
  idx_type ncon = 1;    // number of balance constraints
  idx_type nparts = parts; // number of partitions to create is a parameter
  std::vector<real_type> tpwgts(ncon*nparts, 1./nparts); // load per subdomain and weight (same load on every process)
  std::vector<real_type> ubvec(ncon, 1.05);  // weight tolerance (same weight tolerance for every weight there is)
  idx_type options[3] = {0, 0, 0};                    // use default values for random seed, output and coupling
  idx_type edgecut;                                      // will store number of edges cut by partition
  std::vector<idx_type> part(n);
  
  if (0 == mpihelper.rank()) {
    MPI_Comm comm = Dune::MPIHelper::getLocalCommunicator();
    
#if PARMETIS_MAJOR_VERSION >= 4
    const int OK =
#endif
      
      ParMETIS_V3_PartKway (vtxdist.data(),xadj.data(),adjncy.data(),vwgt.data(),adjwgt.data(),&wgtflag,
			    &numflag,&ncon,&nparts,tpwgts.data(),ubvec.data(),options,&edgecut,part.data(),&comm);  

#if PARMETIS_MAJOR_VERSION >= 4
    if (OK != METIS_OK)
      DUNE_THROW(Dune::Exception, "ParMETIS returned an error code.");
#endif
  }

  std::vector<size_t> retpart(n);
  for (size_t i=0; i<n; i++) retpart[i] = part[i];

  return retpart;
}

/**
 * \brief partition a graph representing subdomains
 *
 * \param subdomain_graph graph represented by std::vector<std::set<size_t>>; a vertex may be coupled to itself
 * \param vertex_weights weight for each vertex in the input graph
 * \param edge_weights weight for each edge in the input graph, only non diagonal entries are used
 * \param mpihelper is required by parmetis
 * \param parts number of partitions in the output
 *
 * Partition the input graph using the method ParMETIS_V3_PartKway of ParMetis. Weights must be given because
 * otherwise the partitioning might be lousy.
 *
 * \return std::vector<unsigned>> for each element stores a set of subomains containing this element
 */
std::vector<size_t> parmetis_sequential_graph_partitioning (const std::vector<std::set<size_t>>& subdomain_graph,
							      const std::vector<size_t> vertex_weights, // for each entry in subdomain_graph
							      const std::vector<std::map<size_t,size_t>> edge_weights, // for each non-diagonal entry in subdomain_graph
							      MPI_Comm comm,
							      int parts)
{
#if PARMETIS_MAJOR_VERSION > 3
  typedef idx_t idx_type;
  typedef ::real_t real_type;
#else
  typedef int idx_type;
  typedef float real_type;
#endif // PARMETIS_MAJOR_VERSION > 3

  // compute some quantities
  size_t nprocs = 1; // number of processors the graph is already distributed on
  size_t n = subdomain_graph.size(); // the number of vertices in the input graph;
  size_t m = 0; // the total number of edges in the input graph
  for (unsigned p=0; p<n; p++)
    for (auto q : subdomain_graph[p])
      if (q!=p)
	m++;

  // prepare ParMetis input
  std::vector<idx_type> vtxdist(nprocs+1); // vertex distribution, everything is sequential here
  vtxdist[0] = 0;
  vtxdist[1] = n;
  std::vector<idx_type> xadj(n+1); // CRS format
  std::vector<idx_type> vwgt(n);
  std::vector<idx_type> adjncy(m);
  std::vector<idx_type> adjwgt(m);
  size_t count=0;
  for (unsigned p=0; p<n; p++)
    {
      xadj[p] = count;
      vwgt[p] = vertex_weights[p];
      for (auto q : subdomain_graph[p])
	if (q!=p)
	  {
	    adjncy[count] = q;
	    adjwgt[count] = edge_weights[p].find(q)->second;
	    count++;
	  }
    }
  xadj[n] = count;
  idx_type wgtflag = 0; // we use weights for vertices and edges
  idx_type numflag = 0; // we are using C-style arrays
  idx_type ncon = 1;    // number of balance constraints
  idx_type nparts = parts; // number of partitions to create is a parameter
  std::vector<real_type> tpwgts(ncon*nparts, 1./nparts); // load per subdomain and weight (same load on every process)
  std::vector<real_type> ubvec(ncon, 1.05);  // weight tolerance (same weight tolerance for every weight there is)
  idx_type options[3] = {0, 0, 0};                    // use default values for random seed, output and coupling
  idx_type edgecut;                                      // will store number of edges cut by partition
  std::vector<idx_type> part(n);

  int rank;
  MPI_Comm_rank (comm, &rank);
  if (rank==0) {
    
#if PARMETIS_MAJOR_VERSION >= 4
    const int OK =
#endif
      
      ParMETIS_V3_PartKway (vtxdist.data(),xadj.data(),adjncy.data(),vwgt.data(),adjwgt.data(),&wgtflag,
			    &numflag,&ncon,&nparts,tpwgts.data(),ubvec.data(),options,&edgecut,part.data(),&comm);  

#if PARMETIS_MAJOR_VERSION >= 4
    if (OK != METIS_OK)
      DUNE_THROW(Dune::Exception, "ParMETIS returned an error code.");
#endif
  }

  std::vector<size_t> retpart(n);
  for (size_t i=0; i<n; i++) retpart[i] = part[i];

  return retpart;
}


#else // HAVE_PARMETIS
#warning "PARMETIS was not found, please check your configuration"
#endif
#endif

/**
 * \brief Takes a partition and extends all subdomains by a number of layers given by overlap
 *
 * \param gv The grid view to be operated on
 * \param partitions number of partitions (or subdomains) the grid has been partitioned into (using the function parmetis_partition)
 * \param partition partition information for each element
 * \param overlap overlap that should be added (zero is fine, then nothing is done)
 * \param mode determines how partitions are grown. Can have the value "vertex" or "element", default is "vertex"
 *
 * If mode has the value "element" the extension is done via element faces. Partitions are extended to neighbors of elements in overlap rounds.
 * If mode has the value "vertex" the extension is done via vertices of the grid. First, partitioning is converted from elements to
 * vertices. Then partitions are extended to neighbors of vertices in overlap rounds. Finally, partitioning
 * is converted back to elements.
 * If mode has neither the value "element" or "vertex" the original partitioning is returned, converted to a set for each element
 *
 * \return std::vector<std::set<unsigned>> for each element stores a set of subomains containing this element
 */
template<class GV>
std::vector<std::set<size_t>> grow_subdomains (const GV& gv, unsigned partitions, const std::vector<size_t>& partition, unsigned overlap, std::string mode="vertex")
{
  const int dim = GV::dimension; // extract dimension (codim of vertices)
  auto& indexset = gv.indexSet(); // to attach data to elements

  if (mode=="vertex")
    {
      std::vector<std::set<size_t>> subdomainsv(indexset.size(dim)); // set of subdomains for each vertex

      // initialize subdomain list for each vertex by the partition
      for (const auto& e : elements(gv))
	for (size_t i=0; i<e.subEntities(dim); ++i)
	  {
	    auto v = e.template subEntity<dim>(i);
	    subdomainsv[indexset.index(v)].insert(partition[indexset.index(e)]);
	  }

      // in each round extend overlap by one
      for (int rounds=0; rounds<overlap; rounds++)
	{
	  std::vector<std::set<size_t>> old(subdomainsv); // copy current state
	  for (const auto& e : elements(gv))
	    {
	      // build union of all partitions in all vertices of the element
	      std::set<size_t> unification;
	      for (size_t i=0; i<e.subEntities(dim); ++i)
		for (const auto& j : old[indexset.index(e.template subEntity<dim>(i))])
		  unification.insert(j);
	      // now add union to all vertices (a clique)
	      for (const auto& j : unification)
		for (size_t i=0; i<e.subEntities(dim); ++i)
		  subdomainsv[indexset.index(e.template subEntity<dim>(i))].insert(j);
	    }
	}

      // now convert again to elements: element is in subdomain if *all* vertices are in subdomain
      std::vector<std::set<size_t>> subdomainse(indexset.size(0)); // set of subdomains for each element
      for (const auto& e : elements(gv))
	{
	  std::set<size_t> intersection(subdomainsv[indexset.index(e.template subEntity<dim>(0))]);
	  for (size_t i=1; i<e.subEntities(dim); ++i)
	    {
	      std::set<size_t> update;
	      for (const auto& j : subdomainsv[indexset.index(e.template subEntity<dim>(i))])
		if (intersection.count(j)>0) update.insert(j);
	      intersection = update;
	    }
	  subdomainse[indexset.index(e)] = intersection;
	}
      // and we are done
      return subdomainse;
    }

  // now the element mode
  std::vector<std::set<size_t>> subdomains(indexset.size(0)); // set of subdomains for each element

  // initialize subdomain list for each element by the partition
  for (const auto& e : elements(gv))
    subdomains[indexset.index(e)].insert(partition[indexset.index(e)]);

  if (mode=="element")
    {
      // in each round extend overlap by one
      for (int rounds=0; rounds<overlap; rounds++)
	{
	  std::vector<std::set<size_t>> old(subdomains); // copy current state
	  for (const auto& e : elements(gv))
	    for (const auto& is : intersections(gv,e))
	      if (is.neighbor())
		for (const auto& i : old[indexset.index(is.outside())])
		  subdomains[indexset.index(e)].insert(i);
	}
    }
  // and we are done
  return subdomains;
}


/**
 * \brief Takes a partition and an input graph and adds overlap layers 
 *
 * \param input_graph graph used to grow overlap
 * \param input_partition partition information for each vertex of the graph
 * \param overlap overlap that should be added (zero is fine, then nothing is done)
 *
 *
 * \return std::vector<std::set<size_t>> for each vertex stores a set of partition numbers containing this vertex
 */
std::vector<std::set<size_t>> sequential_partition_growing (const std::vector<std::set<size_t>>& input_graph,
							    const std::vector<size_t>& input_partition,
							    size_t overlap,
							    const std::vector<bool>& disconnected,
							    bool merge_disconnected)
{
  // size of the input graph 
  auto vertices = input_graph.size();

  // in the output each vertex can be in several partitions
  std::vector<std::set<size_t>> output_partition(vertices); // set of partitions for each vertex

  // in the input each vertex is in exactly one partition
  for (size_t i=0; i<vertices; i++)
    output_partition[i].insert(input_partition[i]);

  // in each round extend overlap by one, zero rounds is fine
  for (int rounds=0; rounds<overlap; rounds++)
    {
      std::vector<std::set<size_t>> old(output_partition); // copy current state
      for (size_t i=0; i<vertices; i++)
	for (auto j : input_graph[i])
	  if (j!=i)
	    for (auto p : old[j])
	      output_partition[i].insert(p);
    }

  if (!merge_disconnected) return output_partition;
  
  // make sure that a disconnected subdomain is not at the boundary
  bool changed=true;
  while (changed)
    {
      changed = false;
      for (size_t i=0; i<vertices; i++) // loop over all vertices of the graph
	for (auto j : input_graph[i]) // loop over all neighbors of this vertex
	  // now j is a neighbor of i
	  for (auto p : output_partition[i])
	    // fine subdomain i is in coarse subdomain p
	    if (output_partition[j].count(p)==0 && disconnected[j])
	      {
		// fine neighboring subdomain j is not in same partition p and it is disconnected; action required
		output_partition[j].insert(p); // add the neighbor also to p
		changed = true; // one more round
		std::cout << "add fine subdomain " << j << " to coarse subdomain " << p << std::endl;
	      }
    }

  return output_partition;
}


#endif // Udune_ftworth_HH
