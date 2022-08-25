#pragma once

#include <cassert>
#include <cstdint>
#include <limits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dune/vtk/types.hh>
#include <dune/vtk/gridcreatorinterface.hh>

namespace Dune
{
  namespace Vtk
  {
    // Create a grid where the input points are not connected and the connectivity
    // describes separated elements.
    template <class Grid>
    struct DiscontinuousGridCreator
        : public GridCreatorInterface<Grid, DiscontinuousGridCreator<Grid>>
    {
      using Self = DiscontinuousGridCreator;
      using Super = GridCreatorInterface<Grid, DiscontinuousGridCreator>;
      using GlobalCoordinate = typename Super::GlobalCoordinate;

      struct CoordLess
      {
        template <class T, int N>
        bool operator() (FieldVector<T,N> const& lhs, FieldVector<T,N> const& rhs) const
        {
          for (int i = 0; i < N; ++i) {
            if (std::abs(lhs[i] - rhs[i]) < std::numeric_limits<T>::epsilon())
              continue;
            return lhs[i] < rhs[i];
          }
          return false;
        }
      };

      using Super::Super;
      using Super::factory;

      void insertVerticesImpl (std::vector<GlobalCoordinate> const& points,
                              std::vector<std::uint64_t> const& /*point_ids*/)
      {
        points_ = &points;
        uniquePoints_.clear();
        std::size_t idx = 0;

        for (auto const& p : points) {
          auto b = uniquePoints_.emplace(std::make_pair(p,idx));
          if (b.second) {
            factory().insertVertex(p);
            ++idx;
          }
        }
      }

      void insertElementsImpl (std::vector<std::uint8_t> const& types,
                              std::vector<std::int64_t> const& offsets,
                              std::vector<std::int64_t> const& connectivity)
      {
        assert(points_ != nullptr);
        std::size_t idx = 0;
        for (std::size_t i = 0; i < types.size(); ++i) {
          auto type = Vtk::to_geometry(types[i]);
          Vtk::CellType cellType{type};

          int nNodes = offsets[i] - (i == 0 ? 0 : offsets[i-1]);
          assert(nNodes > 0);
          std::vector<unsigned int> vtk_cell; vtk_cell.reserve(nNodes);
          for (int j = 0; j < nNodes; ++j) {
            std::size_t v_j = connectivity[idx++];
            std::size_t new_idx = uniquePoints_[(*points_)[v_j]];
            vtk_cell.push_back(new_idx);
          }

          if (cellType.noPermutation()) {
            factory().insertElement(type,vtk_cell);
          } else {
            // apply index permutation
            std::vector<unsigned int> cell(nNodes);
            for (int j = 0; j < nNodes; ++j)
              cell[j] = vtk_cell[cellType.permutation(j)];
            factory().insertElement(type,cell);
          }
        }
      }

    private:
      std::vector<GlobalCoordinate> const* points_ = nullptr;
      std::map<GlobalCoordinate, std::size_t, CoordLess> uniquePoints_;
    };

    // deduction guides
    template <class Grid>
    DiscontinuousGridCreator(GridFactory<Grid>&)
      -> DiscontinuousGridCreator<Grid>;

  } // end namespace Vtk
} // end namespace Dune
