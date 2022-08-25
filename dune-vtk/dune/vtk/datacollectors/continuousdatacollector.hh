#pragma once

#include <numeric>
#include <vector>

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/utility/globalindexset.hh>
#include <dune/vtk/types.hh>

#include "unstructureddatacollector.hh"

namespace Dune
{
  namespace Vtk
  {
    /// Implementation of \ref DataCollector for linear cells, with continuous data.
    template <class GridView, class Partition = Partitions::InteriorBorder>
    class ContinuousDataCollector
        : public UnstructuredDataCollectorInterface<GridView, ContinuousDataCollector<GridView,Partition>, Partition>
    {
      using Self = ContinuousDataCollector;
      using Super = UnstructuredDataCollectorInterface<GridView, Self, Partition>;

    public:
      using Super::dim;
      using Super::partition;

    public:
      ContinuousDataCollector (GridView const& gridView)
        : Super(gridView)
      {}

      /// Collect the vertex indices
      void updateImpl ()
      {
        numPoints_ = 0;
        indexMap_.resize(gridView_.size(dim));
        auto const& indexSet = gridView_.indexSet();
        for (auto const& vertex : vertices(gridView_, partition))
          indexMap_[indexSet.index(vertex)] = std::int64_t(numPoints_++);

        if (gridView_.comm().size() > 1) {
          auto&& e = elements(gridView_, partition);
          numCells_ = std::distance(std::begin(e), std::end(e));
        } else {
          numCells_ = gridView_.size(0);
        }
      }

      /// Return number of grid vertices
      std::uint64_t numPointsImpl () const
      {
        return numPoints_;
      }

      /// Return the coordinates of all grid vertices in the order given by the indexSet
      template <class T>
      std::vector<T> pointsImpl () const
      {
        std::vector<T> data;
        data.reserve(numPoints_ * 3);
        for (auto const& vertex : vertices(gridView_, partition)) {
          auto v = vertex.geometry().center();
          for (std::size_t j = 0; j < v.size(); ++j)
            data.emplace_back(v[j]);
          for (std::size_t j = v.size(); j < 3u; ++j)
            data.emplace_back(0);
        }
        return data;
      }

      /// Return a vector of global unique ids of the points
      std::vector<std::uint64_t> pointIdsImpl () const
      {
        std::vector<std::uint64_t> data;
        data.reserve(numPoints_);
        GlobalIndexSet<GridView> globalIndexSet(gridView_, dim);
        for (auto const& vertex : vertices(gridView_, partition)) {
          data.emplace_back(globalIndexSet.index(vertex));
        }
        return data;
      }

      /// Return number of grid cells
      std::uint64_t numCellsImpl () const
      {
        return numCells_;
      }

      /// Return the types, offsets and connectivity of the cells, using the same connectivity as
      /// given by the grid.
      Cells cellsImpl () const
      {
        auto const& indexSet = gridView_.indexSet();
        auto types = indexSet.types(0);
        int maxVertices = std::accumulate(types.begin(), types.end(), 1, [](int m, GeometryType t) {
          auto refElem = referenceElement<double,dim>(t);
          return std::max(m, refElem.size(dim));
        });

        Cells cells;
        cells.connectivity.reserve(numCells_ * maxVertices);
        cells.offsets.reserve(numCells_);
        cells.types.reserve(numCells_);

        std::int64_t old_o = 0;
        for (auto const& c : elements(gridView_, partition)) {
          Vtk::CellType cellType(c.type());
          for (unsigned int j = 0; j < c.subEntities(dim); ++j)
            cells.connectivity.emplace_back(indexMap_[indexSet.subIndex(c,cellType.permutation(j),dim)]);
          cells.offsets.push_back(old_o += c.subEntities(dim));
          cells.types.push_back(cellType.type());
        }
        return cells;
      }

      /// Evaluate the `fct` at the corners of the elements
      template <class T, class GlobalFunction>
      std::vector<T> pointDataImpl (GlobalFunction const& fct) const
      {
        std::vector<T> data(numPoints_ * fct.numComponents());
        auto const& indexSet = gridView_.indexSet();
        auto localFct = localFunction(fct);
        for (auto const& e : elements(gridView_, partition)) {
          localFct.bind(e);
          Vtk::CellType cellType{e.type()};
          auto refElem = referenceElement(e.geometry());
          for (unsigned int j = 0; j < e.subEntities(dim); ++j) {
            std::size_t idx = fct.numComponents() * indexMap_[indexSet.subIndex(e,cellType.permutation(j),dim)];
            for (int comp = 0; comp < fct.numComponents(); ++comp)
              data[idx + comp] = T(localFct.evaluate(comp, refElem.position(cellType.permutation(j),dim)));
          }
          localFct.unbind();
        }
        return data;
      }

    protected:
      using Super::gridView_;
      std::uint64_t numPoints_ = 0;
      std::uint64_t numCells_ = 0;
      std::vector<std::int64_t> indexMap_;
    };

  } // end namespace Vtk
} // end namespace Dune
