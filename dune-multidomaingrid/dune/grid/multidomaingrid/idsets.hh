#ifndef DUNE_MULTIDOMAINGRID_IDSETS_HH
#define DUNE_MULTIDOMAINGRID_IDSETS_HH

#include <type_traits>

#include <dune/grid/common/indexidset.hh>

namespace Dune {

namespace mdgrid {

template<typename HostGrid, typename MDGridTraits>
class MultiDomainGrid;

template<typename GridImp, typename WrappedIdSet>
class IdSetWrapper :
    public Dune::IdSet<GridImp,IdSetWrapper<GridImp,WrappedIdSet>,
		       typename WrappedIdSet::IdType>
{

  template<typename,typename>
  friend class MultiDomainGrid;

  using HostGrid     = typename std::remove_const_t<GridImp>::HostGrid;
  using Codim0Entity = typename std::remove_const_t<GridImp>::Traits::template Codim<0>::Entity;

public:

  typedef typename WrappedIdSet::IdType IdType;

  template<int codim>
  IdType id(const typename std::remove_const_t<GridImp>::Traits::template Codim<codim>::Entity& e) const {
    return _wrappedIdSet->id(_grid.hostEntity(e));
  }

  template<typename Entity>
  IdType id(const Entity& e) const {
    return _wrappedIdSet->id(_grid.hostEntity(e));
  }

  template<int codim>
  IdType subId(const Codim0Entity& e, int i) const {
    return _wrappedIdSet->subId(_grid.hostEntity(e),i,codim);
  }

  IdType subId(const Codim0Entity& e, int i, unsigned int codim) const {
    return _wrappedIdSet->subId(_grid.hostEntity(e),i,codim);
  }

private:

  const GridImp& _grid;
  const WrappedIdSet* _wrappedIdSet;

  IdSetWrapper(const GridImp& grid) :
    _grid(grid),
    _wrappedIdSet(NULL)
  {}

  void update(const WrappedIdSet& idSet) {
    _wrappedIdSet = &idSet;
  }

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_IDSETS_HH
