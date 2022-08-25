#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <dune/grid/common/gridfactory.hh>
#include <dune/vtk/gridcreatorinterface.hh>
#include <dune/vtk/gridcreators/common.hh>
#include <dune/vtk/gridcreators/continuousgridcreator.hh>

namespace Dune
{
  namespace Vtk
  {
    template <class GridCreator, class Derived>
    struct DerivedGridCreator
        : public GridCreatorInterface<typename GridCreator::Grid, Derived>
    {
      using Self = DerivedGridCreator;
      using Super = GridCreatorInterface<typename GridCreator::Grid, Derived>;
      using Grid = typename GridCreator::Grid;
      using GlobalCoordinate = typename Super::GlobalCoordinate;

      template <class... Args,
        disableCopyMove<DerivedGridCreator, Args...> = 0>
      DerivedGridCreator (Args&&... args)
        : Super(std::forward<Args>(args)...)
        , gridCreator_(Super::factory())
      {}

      void insertVerticesImpl (std::vector<GlobalCoordinate> const& points,
                              std::vector<std::uint64_t> const& point_ids)
      {
        gridCreator_.insertVertices(points, point_ids);
      }

      void insertElementsImpl (std::vector<std::uint8_t> const& types,
                              std::vector<std::int64_t> const& offsets,
                              std::vector<std::int64_t> const& connectivity)
      {
        gridCreator_.insertElements(types, offsets, connectivity);
      }

      void insertPiecesImpl (std::vector<std::string> const& pieces)
      {
        gridCreator_.insertPieces(pieces);
      }

    private:
      GridCreator gridCreator_;
    };

  } // end namespace Vtk;
} // end namespace Dune
