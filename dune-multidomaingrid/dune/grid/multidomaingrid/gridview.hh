#ifndef DUNE_MULTIDOMAINGRID_GRIDVIEW_HH
#define DUNE_MULTIDOMAINGRID_GRIDVIEW_HH

#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/defaultgridview.hh>

namespace Dune {

namespace mdgrid {

template<typename,typename>
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
      typename GridImp::HostGrid::LevelGridView::IntersectionIterator
      >(
        this->grid().hostGrid().levelGridView(entity.level()).ibegin(entity.impl().hostEntity())
        );
  }

  IntersectionIterator iend(const typename BaseT::template Codim<0>::Entity& entity) const
  {
    return IntersectionIteratorWrapper<
      const GridImp,
      typename GridImp::HostGrid::LevelGridView::IntersectionIterator
      >(
        this->grid().hostGrid().levelGridView(entity.level()).iend(entity.impl().hostEntity())
        );
  }

};

template<typename GridImp>
struct LevelGridViewTraits
  : public DefaultLevelGridViewTraits<GridImp>
{
  typedef LevelGridView<GridImp> GridViewImp;
};



template<typename GridImp>
class LeafGridView
  : public DefaultLeafGridView<GridImp>
{

  typedef DefaultLeafGridView<GridImp> BaseT;

public:

  typedef typename BaseT::IntersectionIterator IntersectionIterator;

  LeafGridView(const GridImp& grid)
    : BaseT(grid)
  {}

  IntersectionIterator ibegin(const typename BaseT:: template Codim<0>::Entity& entity) const
  {
    return IntersectionIteratorWrapper<
      const GridImp,
      typename GridImp::HostGrid::LeafGridView::IntersectionIterator
      >(
        this->grid().hostGrid().leafGridView().ibegin(entity.impl().hostEntity())
        );
  }

  IntersectionIterator iend(const typename BaseT:: template Codim<0>::Entity& entity) const
  {
    return IntersectionIteratorWrapper<
      const GridImp,
      typename GridImp::HostGrid::LeafGridView::IntersectionIterator
      >(
        this->grid().hostGrid().leafGridView().iend(entity.impl().hostEntity())
        );
  }

};

template<typename GridImp>
struct LeafGridViewTraits
  : public DefaultLeafGridViewTraits<GridImp>
{
  typedef LeafGridView<GridImp> GridViewImp;
};



} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_GRIDVIEW_HH
