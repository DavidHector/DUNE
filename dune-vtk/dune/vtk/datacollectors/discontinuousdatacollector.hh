#pragma once

#include <vector>

#include <dune/grid/common/partitionset.hh>
#include <dune/vtk/types.hh>

#include "unstructureddatacollector.hh"

namespace Dune
{
  namespace Vtk
  {
    /// Implementation of \ref DataCollector for linear cells, with discontinuous data.
    template <class GridView, class Partition = Partitions::InteriorBorder>
    class DiscontinuousDataCollector
        : public UnstructuredDataCollectorInterface<GridView, DiscontinuousDataCollector<GridView,Partition>, Partition>
    {
      using Self = DiscontinuousDataCollector;
      using Super = UnstructuredDataCollectorInterface<GridView, Self, Partition>;

    public:
      using Super::dim;
      using Super::partition;

    public:
      DiscontinuousDataCollector (GridView const& gridView)
        : Super(gridView)
      {}

      /// Create an index map the uniquely assigns an index to each pair (element,corner)
      void updateImpl ()
      {
        numPoints_ = 0;
        numCells_ = 0;
        indexMap_.resize(gridView_.size(dim));
        std::int64_t vertex_idx = 0;
        auto const& indexSet = gridView_.indexSet();
        for (auto const& c : elements(gridView_, partition)) {
          numCells_++;
          numPoints_ += c.subEntities(dim);
          for (unsigned int i = 0; i < c.subEntities(dim); ++i)
            indexMap_[indexSet.subIndex(c, i, dim)] = vertex_idx++;
        }
      }

      /// The number of points approx. #cell * #corners-per-cell
      std::uint64_t numPointsImpl () const
      {
        return numPoints_;
      }

      /// Return the coordinates of the corners of all cells
      template <class T>
      std::vector<T> pointsImpl () const
      {
        std::vector<T> data(numPoints_ * 3);
        auto const& indexSet = gridView_.indexSet();
        for (auto const& element : elements(gridView_, partition)) {
          for (unsigned int i = 0; i < element.subEntities(dim); ++i) {
            std::size_t idx = 3 * indexMap_[indexSet.subIndex(element, i, dim)];
            auto v = element.geometry().corner(i);
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
        return numCells_;
      }

      /// Connect the corners of each cell. The leads to a global discontinuous grid
      Cells cellsImpl () const
      {
        Cells cells;
        cells.connectivity.reserve(numPoints_);
        cells.offsets.reserve(numCells_);
        cells.types.reserve(numCells_);

        std::int64_t old_o = 0;
        auto const& indexSet = gridView_.indexSet();
        for (auto const& c : elements(gridView_, partition)) {
          Vtk::CellType cellType(c.type());
          for (unsigned int j = 0; j < c.subEntities(dim); ++j) {
            std::int64_t vertex_idx = indexMap_[indexSet.subIndex(c,cellType.permutation(j),dim)];
            cells.connectivity.push_back(vertex_idx);
          }
          cells.offsets.push_back(old_o += c.subEntities(dim));
          cells.types.push_back(cellType.type());
        }

        return cells;
      }

      /// Evaluate the `fct` in the corners of each cell
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
            std::size_t idx = fct.numComponents() * indexMap_[indexSet.subIndex(e, cellType.permutation(j), dim)];
            for (int comp = 0; comp < fct.numComponents(); ++comp)
              data[idx + comp] = T(localFct.evaluate(comp, refElem.position(cellType.permutation(j),dim)));
          }
          localFct.unbind();
        }
        return data;
      }

    protected:
      using Super::gridView_;
      std::uint64_t numCells_ = 0;
      std::uint64_t numPoints_ = 0;
      std::vector<std::int64_t> indexMap_;
    };

  } // end namespace Vtk
} // end namespace Dune
