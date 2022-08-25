#pragma once

#include <cassert>
#include <map>
#include <vector>

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/vtk/types.hh>
#include <dune/vtk/utility/lagrangepoints.hh>

#include "unstructureddatacollector.hh"

namespace Dune
{
  namespace Vtk
  {
    /// Implementation of \ref DataCollector for lagrange cells
    template <class GridView, int ORDER = -1>
    class LagrangeDataCollector
        : public UnstructuredDataCollectorInterface<GridView, LagrangeDataCollector<GridView,ORDER>, Partitions::All>
    {
      using Self = LagrangeDataCollector;
      using Super = UnstructuredDataCollectorInterface<GridView, Self, Partitions::All>;

    public:
      static_assert(ORDER != 0, "Order 0 not supported");
      using Super::dim;
      using Super::partition; // NOTE: lagrange data-collector currently implemented for the All partition only

    public:
      LagrangeDataCollector (GridView const& gridView, int order = ORDER)
        : Super(gridView)
        , order_(order)
      {
        assert(order > 0 && "Order 0 not supported");
        assert(ORDER < 0 || order == ORDER);
      }

      /// Construct the point sets
      void updateImpl ()
      {
        auto const& indexSet = gridView_.indexSet();

        pointSets_.clear();
        for (auto gt : indexSet.types(0))
          pointSets_.emplace(gt, order_);

        for (auto& pointSet : pointSets_)
          pointSet.second.build(pointSet.first);

        numPoints_ = indexSet.size(dim);
        for (auto const& pointSet : pointSets_) {
          auto gt = pointSet.first;
          auto refElem = referenceElement<double,dim>(gt);
          numPoints_ += (pointSet.second.size() - refElem.size(dim)) * indexSet.size(gt);
        }
      }

      /// Return number of lagrange nodes
      std::uint64_t numPointsImpl () const
      {
        return numPoints_;
      }

      /// Return a vector of point coordinates.
      /**
      * The vector of point coordinates is composed of vertex coordinates first and second
      * edge center coordinates.
      **/
      template <class T>
      std::vector<T> pointsImpl () const
      {
        std::vector<T> data(this->numPoints() * 3);
        auto const& indexSet = gridView_.indexSet();

        std::size_t shift = indexSet.size(dim);

        for (auto const& element : elements(gridView_, partition)) {
          auto geometry = element.geometry();
          auto refElem = referenceElement<T,dim>(element.type());

          auto const& pointSet = pointSets_.at(element.type());
          unsigned int vertexDOFs = refElem.size(dim);
          unsigned int innerDOFs = pointSet.size() - vertexDOFs;

          for (std::size_t i = 0; i < pointSet.size(); ++i) {
            auto const& p = pointSet[i];
            if (i < vertexDOFs)
              assert(p.localKey().codim() == dim);

            auto const& localKey = p.localKey();
            std::size_t idx = 3 * (localKey.codim() == dim
                ? indexSet.subIndex(element, localKey.subEntity(), dim)
                : innerDOFs*indexSet.index(element) + (i - vertexDOFs) + shift);

            auto v = geometry.global(p.point());
            for (std::size_t j = 0; j < v.size(); ++j)
              data[idx + j] = T(v[j]);
            for (std::size_t j = v.size(); j < 3u; ++j)
              data[idx + j] = T(0);
          }
        }
        return data;
      }

      /// Return number of grid cells
      std::uint64_t numCellsImpl () const
      {
        return gridView_.size(0);
      }

      /// \brief Return cell types, offsets, and connectivity. \see Cells
      /**
      * The cell connectivity is composed of cell vertices first and second cell edges,
      * where the indices are grouped [vertex-indices..., (#vertices)+edge-indices...]
      **/
      Cells cellsImpl () const
      {
        Cells cells;
        cells.connectivity.reserve(this->numPoints());
        cells.offsets.reserve(this->numCells());
        cells.types.reserve(this->numCells());

        auto const& indexSet = gridView_.indexSet();
        std::size_t shift = indexSet.size(dim);

        std::int64_t old_o = 0;
        for (auto const& element : elements(gridView_, partition)) {
          auto refElem = referenceElement<double,dim>(element.type());
          Vtk::CellType cellType(element.type(), Vtk::CellType::LAGRANGE);

          auto const& pointSet = pointSets_.at(element.type());
          unsigned int vertexDOFs = refElem.size(dim);
          unsigned int innerDOFs = pointSet.size() - vertexDOFs;

          for (std::size_t i = 0; i < pointSet.size(); ++i) {
            auto const& p = pointSet[i];
            auto const& localKey = p.localKey();
            std::size_t idx = (localKey.codim() == dim
                ? indexSet.subIndex(element, localKey.subEntity(), dim)
                : innerDOFs*indexSet.index(element) + (i - vertexDOFs) + shift);
            cells.connectivity.push_back(idx);
          }

          cells.offsets.push_back(old_o += pointSet.size());
          cells.types.push_back(cellType.type());
        }
        return cells;
      }

      /// Evaluate the `fct` at element vertices and edge centers in the same order as the point coords.
      template <class T, class GlobalFunction>
      std::vector<T> pointDataImpl (GlobalFunction const& fct) const
      {
        int nComps = fct.numComponents();
        std::vector<T> data(this->numPoints() * nComps);
        auto const& indexSet = gridView_.indexSet();

        std::size_t shift = indexSet.size(dim);

        auto localFct = localFunction(fct);
        for (auto const& element : elements(gridView_, partition)) {
          localFct.bind(element);
          auto refElem = referenceElement<T,dim>(element.type());

          auto const& pointSet = pointSets_.at(element.type());
          unsigned int vertexDOFs = refElem.size(dim);
          unsigned int innerDOFs = pointSet.size() - vertexDOFs;

          for (std::size_t i = 0; i < pointSet.size(); ++i) {
            auto const& p = pointSet[i];
            auto const& localKey = p.localKey();
            std::size_t idx = nComps * (localKey.codim() == dim
                ? indexSet.subIndex(element, localKey.subEntity(), dim)
                : innerDOFs*indexSet.index(element) + (i - vertexDOFs) + shift);

            for (int comp = 0; comp < nComps; ++comp)
              data[idx + comp] = T(localFct.evaluate(comp, p.point()));
          }
          localFct.unbind();
        }
        return data;
      }

    protected:
      using Super::gridView_;

      unsigned int order_;
      std::uint64_t numPoints_ = 0;

      using PointSet = LagrangePointSet<typename GridView::ctype, GridView::dimension>;
      std::map<GeometryType, PointSet> pointSets_;
    };

  } // end namespace Vtk
} // end namespace Dune
