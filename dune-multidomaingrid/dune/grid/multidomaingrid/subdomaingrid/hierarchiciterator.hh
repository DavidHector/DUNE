#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_HIERARCHICITERATOR_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_HIERARCHICITERATOR_HH

#include "entity.hh"

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<int,int,typename>
class EntityWrapper;

template<typename GridImp>
class HierarchicIteratorWrapper
{

  template< int cd, class Grid, class IteratorImp >
  friend class Dune::EntityIterator;

  template<int, int, typename>
  friend class EntityWrapper;

  template<int, int, typename>
  friend class EntityWrapper;

  using MultiDomainIterator = typename GridImp::MultiDomainGrid::template Codim<0>::Entity::HierarchicIterator;

public:

  static const int codimension = 0;

  using EntityWrapper        = Dune::mdgrid::subdomain::EntityWrapper<0,GridImp::dimension,GridImp>;
  using Entity               = typename GridImp::template Codim<0>::Entity;

  HierarchicIteratorWrapper()
    : _grid(nullptr)
  {}

  HierarchicIteratorWrapper(
    const GridImp* grid,
    const MultiDomainIterator& multiDomainIterator,
    const MultiDomainIterator& multiDomainEnd
    )
    : _grid(grid)
    , _multiDomainIterator(multiDomainIterator)
    , _multiDomainEnd(multiDomainEnd)
  {
    incrementToNextValidPosition();
  }

  bool equals(const HierarchicIteratorWrapper& r) const
  {
    return _grid == r._grid && _multiDomainIterator == r._multiDomainIterator;
  }

  Entity dereference() const
  {
    return {EntityWrapper(_grid,*_multiDomainIterator)};
  }

  void incrementToNextValidPosition() {
    while(_multiDomainIterator != _multiDomainEnd && !_grid->containsMultiDomainEntity(*_multiDomainIterator))
      {
        ++_multiDomainIterator;
      }
  }

  void increment() {
    ++_multiDomainIterator;
    incrementToNextValidPosition();
  }

private:

  const GridImp* _grid;
  MultiDomainIterator _multiDomainIterator;
  MultiDomainIterator _multiDomainEnd;

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_HIERARCHICITERATOR_HH
