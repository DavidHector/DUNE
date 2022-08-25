#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <dune/grid/common/gridfactory.hh>
#include <dune/vtk/gridcreatorinterface.hh>
#include <dune/vtk/gridcreators/common.hh>
#include <dune/vtk/gridcreators/derivedgridcreator.hh>
#include <dune/vtk/gridcreators/continuousgridcreator.hh>

namespace Dune
{
  namespace Vtk
  {
    // create a distributed grid in parallel. Currently only supported by ALUGrid
    template <class Grid>
    struct ParallelGridCreator
        : public DerivedGridCreator<ContinuousGridCreator<Grid>, ParallelGridCreator<Grid>>
    {
      using Self = ParallelGridCreator;
      using Super = DerivedGridCreator<ContinuousGridCreator<Grid>, Self>;
      using GlobalCoordinate = typename Super::GlobalCoordinate;
      using VertexId = VertexId_t<GridFactory<Grid>>;

      // The GridFactory must support insertion of global vertex IDs
      static_assert(Std::is_detected<HasInsertVertex, GridFactory<Grid>, GlobalCoordinate, VertexId>{}, "");

    public:

      using Super::Super;

      void insertVerticesImpl (std::vector<GlobalCoordinate> const& points,
                               std::vector<std::uint64_t> const& point_ids)
      {
        assert(point_ids.size() == points.size());
        for (std::size_t i = 0; i < points.size(); ++i)
          this->factory().insertVertex(points[i], VertexId(point_ids[i]));
      }

      void insertPiecesImpl (std::vector<std::string> const& pieces)
      {
        if (int(pieces.size()) == this->comm().size()) {
          VtkReader<Grid, Self> pieceReader(this->factory());
          pieceReader.read(pieces[this->comm().rank()], true);
        }
      }
    };

    // deduction guides
    template <class Grid>
    ParallelGridCreator(GridFactory<Grid>&)
      -> ParallelGridCreator<Grid>;

  } // end namespace Vtk
} // end namespace Dune
