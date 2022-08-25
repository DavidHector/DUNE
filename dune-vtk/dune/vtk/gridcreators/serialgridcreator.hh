#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <dune/grid/common/gridfactory.hh>
#include <dune/vtk/gridcreatorinterface.hh>
#include <dune/vtk/gridcreators/discontinuousgridcreator.hh>

namespace Dune
{
  namespace Vtk
  {
    // create a distributed grid on rank 0. Needs to be load balanced afterwards.
    template <class Grid>
    struct SerialGridCreator
        : public GridCreatorInterface<Grid, SerialGridCreator<Grid>>
    {
      using Self = SerialGridCreator;
      using Super = GridCreatorInterface<Grid, Self>;
      using GlobalCoordinate = typename Super::GlobalCoordinate;

    public:

      using Super::Super;

      void insertVerticesImpl (std::vector<GlobalCoordinate> const& points,
                               std::vector<std::uint64_t> const& /*point_ids*/)
      {
        shift_.push_back(points_.size());
        points_.reserve(points_.size() + points.size());
        points_.insert(points_.end(), points.begin(), points.end());
      }

      void insertElementsImpl (std::vector<std::uint8_t> const& types,
                               std::vector<std::int64_t> const& offsets,
                               std::vector<std::int64_t> const& connectivity)
      {
        types_.reserve(types_.size() + types.size());
        types_.insert(types_.end(), types.begin(), types.end());

        offsets_.reserve(offsets_.size() + offsets.size());
        std::transform(offsets.begin(), offsets.end(), std::back_inserter(offsets_),
          [shift=offsets_.empty() ? 0 : offsets_.back()](std::int64_t o) { return o + shift; });

        connectivity_.reserve(connectivity_.size() + connectivity.size());
        std::transform(connectivity.begin(), connectivity.end(), std::back_inserter(connectivity_),
          [shift=shift_.back()](std::int64_t idx) { return idx + shift; });
      }

      void insertPiecesImpl (std::vector<std::string> const& pieces)
      {
        if (this->comm().rank() == 0) {
          VtkReader<Grid, Self> pieceReader(*this);
          for (std::string const& piece : pieces) {
            pieceReader.read(piece, false);
            pieceReader.fillGridCreator(false);
          }

          DiscontinuousGridCreator<Grid> creator(this->factory());
          creator.insertVertices(points_, {});
          creator.insertElements(types_, offsets_, connectivity_);
        }
      }

    private:
      std::vector<GlobalCoordinate> points_;
      std::vector<std::uint8_t> types_;
      std::vector<std::int64_t> offsets_;
      std::vector<std::int64_t> connectivity_;
      std::vector<std::int64_t> shift_;
    };

    // deduction guides
    template <class Grid>
    SerialGridCreator(GridFactory<Grid>&)
      -> SerialGridCreator<Grid>;

  } // end namespace Vtk
} // end namespace Dune
