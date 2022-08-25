#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INDEXSETS_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INDEXSETS_HH

#include <unordered_map>
#include <vector>
#include <array>
#include <algorithm>
#include <type_traits>
#include <tuple>

#include <dune/grid/common/indexidset.hh>
#include <dune/grid/multidomaingrid/utility.hh>

#include <type_traits>

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<typename MDGrid>
class SubDomainGrid;


template<typename GridImp, typename MDIndexSet>
class IndexSetWrapper :
    public Dune::IndexSet<GridImp,IndexSetWrapper<GridImp,MDIndexSet>,
                          typename MDIndexSet::IndexType,
                          typename MDIndexSet::Types>
{

  template<typename, typename>
  friend class IndexSetWrapper;

  template<typename MDGrid>
  friend class SubDomainGrid;

  template<typename, typename, typename>
  friend class IntersectionWrapper;

  template<typename, typename, int, PartitionIteratorType, typename>
  friend class IteratorWrapper;

  typedef IndexSetWrapper<GridImp,MDIndexSet> ThisType;

  using HostGrid = typename std::remove_const_t<GridImp>::HostGrid;
  using MDGrid   = typename std::remove_const_t<GridImp>::MultiDomainGrid;
  using ctype    = typename std::remove_const_t<GridImp>::ctype;

public:

  typedef typename MDIndexSet::Types Types;

  using SubDomainIndex = typename std::remove_const_t<GridImp>::SubDomainIndex;
  typedef typename MDIndexSet::IndexType IndexType;
  static const int dimension = std::remove_const_t<GridImp>::dimension;

private:

  using Codim0Entity = typename std::remove_const_t<GridImp>::Traits::template Codim<0>::Entity;

public:

  template<int codim>
  IndexType index(const typename std::remove_const_t<GridImp>::Traits::template Codim<codim>::Entity& e) const {
    return _mdIndexSet.template index<codim>(_grid.domain(),_grid.multiDomainEntity(e));
  }

  template<typename Entity>
  IndexType index(const Entity& e) const {
    return index<Entity::codimension>(e);
  }

  template<int codim, typename Entity>
  IndexType subIndex(const Entity& e, int i) const {
    return _mdIndexSet.subIndex(_grid.domain(),_grid.multiDomainEntity(e),i,codim);
  }

  template<typename Entity>
  IndexType subIndex(const Entity& e, int i, unsigned int codim) const {
    return _mdIndexSet.subIndex(_grid.domain(),_grid.multiDomainEntity(e),i,codim);
  }

  Types types(int codim) const {
    return _mdIndexSet.types(_grid.domain(),codim);
  }

  IndexType size(GeometryType type) const {
    return _mdIndexSet.sizeForSubDomain(_grid.domain(),type);
  }

  IndexType size(int codim) const {
    return _mdIndexSet.sizeForSubDomain(_grid.domain(),codim);
  }

  template<typename EntityType>
  bool contains(const EntityType& e) const {
    return _mdIndexSet.contains(_grid.domain(),_grid.multiDomainEntity(e));
  }

  bool operator==(const IndexSetWrapper& rhs) const {
    return (_grid == rhs._grid && &_mdIndexSet == &rhs._mdIndexSet);
  }

  IndexSetWrapper(const GridImp& grid, const MDIndexSet& mdIndexSet) :
    _grid(grid),
    _mdIndexSet(mdIndexSet)
  {}

private:

  const GridImp& _grid;
  const MDIndexSet& _mdIndexSet;

  template<typename EntityType>
  bool containsHostEntity(const EntityType& he) const {
    return _mdIndexSet.containsForSubDomain(_grid.domain(),he);
  }

  template<typename EntityType>
  bool containsMultiDomainEntity(const EntityType& mde) const {
    return _mdIndexSet.contains(_grid.domain(),mde);
  }

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INDEXSETS_HH
