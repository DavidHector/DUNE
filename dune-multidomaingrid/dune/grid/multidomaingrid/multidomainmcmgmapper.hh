#ifndef DUNE_MULTIDOMAINGRID_MULTDIDOMAINMCMGMAPPER_HH
#define DUNE_MULTIDOMAINGRID_MULTDIDOMAINMCMGMAPPER_HH

#include <iostream>
#include <map>
#include <numeric>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/referenceelements.hh>

namespace Dune {

namespace mdgrid {

/**
 * @addtogroup Mapper
 *
 * @{
 */

namespace {

  template<typename GV, bool max_subomain_index_is_static>
  struct MCMGMapperStorageProvider
  {

  protected:

    typedef typename GV::IndexSet::IndexType IndexType;

    void update(const GV& gv) {}

    std::array<std::vector<IndexType>,GV::Grid::maxSubDomainIndex()> _offsets; // provide a map with all geometry types

  };

  template<typename GV>
  struct MCMGMapperStorageProvider<GV,false>
  {

  protected:

    typedef typename GV::IndexSet::IndexType IndexType;

    void update(const GV& gv) {
      _offsets.resize(gv.grid().maxSubDomainIndex());
    }

    std::vector<std::vector<IndexType> > _offsets; // provide a map with all geometry types

  };


}

/** @brief Implementation class for a multiple codim and multiple geometry type mapper.
 *
 * In this implementation of a mapper the entity set used as domain for the map consists
 * of the entities of a subset of codimensions in the given index set. The index
 * set may contain entities of several geometry types. This
 * version is usually not used directly but is used to implement versions for leafwise and levelwise
 * entity sets.
 *
 * Template parameters are:
 *
 * \par GV
 *    A Dune GridView type.
 * \par Layout
 *  A helper class with a method contains(), that returns true for all geometry
 *  types that are in the domain of the map.  The class should be of the following
 *  shape
 \code
 template<int dim>
 struct LayoutClass {
 bool contains (Dune::GeometryType gt) const {
 // Return true if gt is in the domain of the map
 }
 };
 \endcode
 *
 * If you don't want to use the default constructor of the LayoutClass you can construct it yourself
 * and hand it to the respective constructor.
 */
template <typename GV>
class MultiDomainMCMGMapper : public MultipleCodimMultipleGeomTypeMapper<GV>,
                              public MCMGMapperStorageProvider<GV,GV::Grid::maxSubDomainIndexIsStatic()>

{

  typedef MultipleCodimMultipleGeomTypeMapper<GV> Base;

  typedef MCMGMapperStorageProvider<GV,GV::Grid::maxSubDomainIndexIsStatic()> StorageProvider;

  using StorageProvider::_offsets;

public:

  typedef typename GV::IndexSet::IndexType IndexType;
  typedef typename GV::Grid::SubDomainIndex SubDomainIndex;

  MultiDomainMCMGMapper (const GV& gv, const MCMGLayout& layout)
    : Base(gv,layout)
  {
    update(gv);
  }

  /** @brief Construct mapper from grid and one of its index sets.

      \param grid A Dune grid object.
      \param indexset IndexSet object returned by grid.

  */
  MultiDomainMCMGMapper (const GV& gv)
    : Base(gv)
  {
    update(gv);
  }

  using Base::gridView;

  /** @brief Map entity to array index.

      \param e Reference to codim cc entity, where cc is the template parameter of the function.
      \return An index in the range 0 ... Max number of entities in set - 1.
  */
  template<class EntityType>
  int map (SubDomainIndex subDomain, const EntityType& e) const
  {
    return gridView().indexSet().index(subDomain,e) + _offsets[subDomain][GlobalGeometryTypeIndex::index(e.type())];
  }

  /** @brief Map subentity of codim 0 entity to array index.

      \param e Reference to codim 0 entity.
      \param i Number of subentity of e
      \param codim Codimension of the subendity
      \return An index in the range 0 ... Max number of entities in set - 1.
  */
  int map (SubDomainIndex subDomain, const typename GV::template Codim<0>::Entity& e, int i, unsigned int codim) const
  {
    GeometryType gt =
        ReferenceElements<double, GV::dimension>::general(e.type()).type(i,
                                                                         codim);
    return gridView().indexSet().subIndex(subDomain, e, i, codim) +
           _offsets[subDomain][GlobalGeometryTypeIndex::index(gt)];
  }

  /** @brief Return total number of entities in the entity set managed by the mapper.

      This number can be used to allocate a vector of data elements associated with the
      entities of the set. In the parallel case this number is per process (i.e. it
      may be different in different processes).

      \return Size of the entity set.
  */
  int size (SubDomainIndex subDomain) const
  {
    return _offsets[subDomain].back();
  }

  /** @brief Returns true if the entity is contained in the index set

      \param e Reference to entity
      \param result integer reference where corresponding index is  stored if true
      \return true if entity is in entity set of the mapper
  */
  template<class EntityType>
  bool contains (SubDomainIndex subDomain, const EntityType& e, IndexType& result) const
  {
    if (!gridView().indexSet().contains(subDomain, e) ||
        !Base::layout()(e.type(), GV::dimension)) {
      result = 0;
      return false;
    }
    result = map(subDomain,e);
    return true;
  }

  /** @brief Returns true if the entity is contained in the index set

      \param e Reference to codim 0 entity
      \param i subentity number
      \param result integer reference where corresponding index is  stored if true
      \return true if entity is in entity set of the mapper
  */
  template<int cc> // this is now the subentity's codim
  bool contains (SubDomainIndex subDomain, const typename GV::template Codim<0>::Entity& e, int i, IndexType& result) const
  {
    GeometryType gt = ReferenceElements<double,GV::dimension>::general(e.type()).type(i,cc);
    // if the entity is contained in the subdomain, all of its subentities are contained as well
    if (!gridView().indexSet().contains(subDomain, e) ||
        !Base::layout()(gt, GV::dimension)) {
      result = 0;
      return false;
    }
    result = gridView().indexSet().subIndex(subDomain, e, i, cc) +
             _offsets[subDomain][GlobalGeometryTypeIndex::index(gt)];
    return true;
  }

  /** @brief Recalculates indices after grid adaptation
   */
  [[deprecated("Use update(gridView) instead! Will be removed after release 2.8.")]]
  void update ()
  {
    update(gridView());
  }

  /** @brief Recalculates indices after grid adaptation
   *
   * After grid adaptation you need to call this to update
   * the stored gridview and recalculate the indices.
   */
  void update(const GV& gv)
  {
    static_cast<Base*>(this)->update(gv);
    StorageProvider::update(gv);
    for (SubDomainIndex subDomain = 0;
         subDomain < gridView().grid().maxSubDomainIndex(); ++subDomain) {
      std::vector<IndexType>& offsets = _offsets[subDomain];
      offsets.resize(GlobalGeometryTypeIndex::size(GV::dimension) + 1);
      std::fill(offsets.begin(),offsets.end(),0);

      // Compute offsets for the different geometry types.
      // Note that mapper becomes invalid when the grid is modified.
      for (int cc = 0; cc <= GV::dimension; ++cc)
        {
          for (auto gt : gridView().indexSet().types(subDomain,cc))
            offsets[GlobalGeometryTypeIndex::index(gt) + 1] =
                gridView().indexSet().size(subDomain, gt);
          // convert sizes to offset
          // last entry stores total size
          std::partial_sum(offsets.begin(),offsets.end(),offsets.begin());
        }
    }
  }
};

/** @} */

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_MULTDIDOMAINMCMGMAPPER_HH
