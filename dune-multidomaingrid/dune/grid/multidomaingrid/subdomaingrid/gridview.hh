#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_GRIDVIEW_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_GRIDVIEW_HH

#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/defaultgridview.hh>

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<typename,typename,typename>
class IntersectionIteratorWrapper;

template<typename GridImp>
class LevelGridView
  : public DefaultLevelGridView<GridImp>
{

  using BaseT = DefaultLevelGridView<GridImp>;

public:

  using typename BaseT::IntersectionIterator;

  LevelGridView(const GridImp& grid, int level)
    : BaseT(grid,level)
  {}

  IntersectionIterator ibegin(const typename BaseT::template Codim<0>::Entity& entity) const
  {
    return IntersectionIteratorWrapper<
      const GridImp,
      typename BaseT::IndexSet,
      typename GridImp::MultiDomainGrid::LevelGridView::IntersectionIterator
      >(
        &this->indexSet(),
        this->grid()._grid.levelGridView(entity.level()).ibegin(entity.impl().multiDomainEntity())
        );
  }

  IntersectionIterator iend(const typename BaseT::template Codim<0>::Entity& entity) const
  {
    return IntersectionIteratorWrapper<
      const GridImp,
      typename BaseT::IndexSet,
      typename GridImp::MultiDomainGrid::LevelGridView::IntersectionIterator
      >(
        &this->indexSet(),
        this->grid()._grid.levelGridView(entity.level()).iend(entity.impl().multiDomainEntity())
        );
  }

};

template<typename GridImp>
struct LevelGridViewTraits
  : public DefaultLevelGridViewTraits<GridImp>
{
  using GridViewImp = LevelGridView<GridImp>;
};



template<typename GridImp>
class LeafGridView
  : public DefaultLeafGridView<GridImp>
{

  using BaseT = DefaultLeafGridView<GridImp>;

public:

  using typename BaseT::IntersectionIterator;

  LeafGridView(const GridImp& grid)
    : BaseT(grid)
  {}

  IntersectionIterator ibegin(const typename BaseT::template Codim<0>::Entity& entity) const
  {
    return IntersectionIteratorWrapper<
      const GridImp,
      typename BaseT::IndexSet,
      typename GridImp::MultiDomainGrid::LeafGridView::IntersectionIterator
      >(
        &this->indexSet(),
        this->grid()._grid.leafGridView().ibegin(entity.impl().multiDomainEntity())
        );
  }

  IntersectionIterator iend(const typename BaseT::template Codim<0>::Entity& entity) const
  {
    return IntersectionIteratorWrapper<
      const GridImp,
      typename BaseT::IndexSet,
      typename GridImp::MultiDomainGrid::LeafGridView::IntersectionIterator
      >(
        &this->indexSet(),
        this->grid()._grid.leafGridView().iend(entity.impl().multiDomainEntity())
        );
  }

};

template<typename GridImp>
struct LeafGridViewTraits
  : public DefaultLeafGridViewTraits<GridImp>
{
  using GridViewImp = LeafGridView<GridImp>;
};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_GRIDVIEW_HH
