#ifndef DUNE_GRID_MULTIDOMAINGRID_SUBDOMAINGRID_ITERATOR_HH
#define DUNE_GRID_MULTIDOMAINGRID_SUBDOMAINGRID_ITERATOR_HH

#include <dune/grid/common/gridenums.hh>

#include "entity.hh"

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<typename MDGrid>
class SubDomainGrid;

template<int,int,typename>
class EntityWrapper;

template<
  typename GridView,
  typename MultiDomainIterator,
  int codim,
  PartitionIteratorType pitype,
  typename GridImp
  >
class IteratorWrapper
{

  template<typename MDGrid>
  friend class SubDomainGrid;

  template< int cd, class Grid, class IteratorImp >
  friend class Dune::EntityIterator;

  static const int codimension = codim;

  using IndexSet      = typename GridView::IndexSet;

public:

  // Entity has to be public for iterator traits
  using Entity        = typename GridImp::template Codim<codim>::Entity;

private:

  using EntityWrapper        = Dune::mdgrid::subdomain::EntityWrapper<codim,GridImp::dimension,GridImp>;

  IteratorWrapper()
    : _grid(nullptr)
    , _indexSet(nullptr)
  {}

  IteratorWrapper(
    const GridImp* grid,
    const IndexSet* indexSet,
    MultiDomainIterator multiDomainIterator,
    MultiDomainIterator endIterator
    )
    : _grid(grid)
    , _indexSet(indexSet)
    , _multiDomainIterator(multiDomainIterator)
    , _end(endIterator)
  {
    incrementToNextValidPosition();
  }

  void incrementToNextValidPosition() {
    while(_multiDomainIterator != _end && !_indexSet->containsMultiDomainEntity(*_multiDomainIterator))
      {
        ++_multiDomainIterator;
      }
  }

  void increment() {
    ++_multiDomainIterator;
    incrementToNextValidPosition();
  }

  bool equals(const IteratorWrapper& r) const
  {
    return _grid == r._grid && _indexSet == r._indexSet && _multiDomainIterator == r._multiDomainIterator;
  }

  Entity dereference() const
  {
    return {EntityWrapper(_grid,*_multiDomainIterator)};
  }

  int level() const
  {
    return _multiDomainIterator.level();
  }

public:

private:

  const GridImp* _grid;
  const IndexSet* _indexSet;
  MultiDomainIterator _multiDomainIterator;
  MultiDomainIterator _end;

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_GRID_MULTIDOMAINGRID_SUBDOMAINGRID_ITERATOR_HH
