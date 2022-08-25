#pragma once

#include <array>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/vtk/datacollectorinterface.hh>

#include "continuousdatacollector.hh"

namespace Dune
{
  namespace Vtk
  {
    /// The Interface for structured data-collectors
    template <class GridView, class Derived>
    class StructuredDataCollectorInterface
        : public DataCollectorInterface<GridView, Derived>
    {
    protected:
      using Super = DataCollectorInterface<GridView, Derived>;
      using SubDataCollector = ContinuousDataCollector<GridView>;
      using ctype = typename GridView::ctype;

    public:
      using Super::dim;
      using Super::partition;

    public:
      StructuredDataCollectorInterface (GridView const& gridView)
        : Super(gridView)
        , subDataCollector_(gridView)
      {}

      /// Sequence of Index pairs [begin, end) for the cells in each direction
      std::array<int, 6> wholeExtent () const
      {
        return this->asDerived().wholeExtentImpl();
      }

      /// Sequence of Index pairs [begin, end) for the cells in each direction of the local partition
      std::array<int, 6> extent () const
      {
        return this->asDerived().extentImpl();
      }

      /// Call the `writer` with extent
      template <class Writer>
      void writeLocalPiece (Writer const& writer) const
      {
        this->asDerived().writeLocalPieceImpl(writer);
      }

      /// Call the `writer` with piece number and piece extent
      template <class Writer>
      void writePieces (Writer const& writer) const
      {
        this->asDerived().writePiecesImpl(writer);
      }

      /// Interface for ImageData:
      /// @{

      /// Lower left corner of the grid
      FieldVector<ctype, 3> origin () const
      {
        return this->asDerived().originImpl();
      }

      /// Constant grid spacing in each coordinate direction
      FieldVector<ctype, 3> spacing () const
      {
        return this->asDerived().spacingImpl();
      }

      /// @}

      /// Interface for RectilinearGrid
      /// @{

      /// The coordinates defines point coordinates for an extent by specifying the ordinate along each axis.
      template <class T>
      std::array<std::vector<T>, 3> coordinates () const
      {
        return this->asDerived().template coordinatesImpl<T>();
      }

      /// @}


    public: // default implementation:

      /// \copyref DefaultDataCollector::update.
      void updateImpl ()
      {
        subDataCollector_.update();

    #if HAVE_MPI
        int rank = gridView_.comm().rank();
        int numRanks = gridView_.comm().size();

        if (rank == 0) {
          extents_.resize(numRanks);
          requests_.resize(numRanks, MPI_REQUEST_NULL);
          for (int i = 1; i < numRanks; ++i)
            MPI_Irecv(extents_[i].data(), extents_[i].size(), MPI_INT, i, /*tag=*/6, gridView_.comm(), &requests_[i]);
        }

        sendRequest_ = MPI_REQUEST_NULL;
    #endif
      }

      /// Return number of grid cells
      std::uint64_t numCellsImpl () const
      {
        auto extent = this->extent();
        std::uint64_t num = 1;
        for (int d = 0; d < dim; ++d)
          num *= extent[2*d+1] - extent[2*d];
        return num;
      }

      /// Return number of grid vertices
      std::uint64_t numPointsImpl () const
      {
        return subDataCollector_.numPoints();
      }

      /// \copydoc DefaultDataCollector::points.
      template <class T>
      std::vector<T> pointsImpl () const
      {
        return subDataCollector_.template points<T>();
      }

      /// \copydoc DefaultDataCollector::pointData
      template <class T, class GlobalFunction>
      std::vector<T> pointDataImpl (GlobalFunction const& fct) const
      {
        return subDataCollector_.template pointData<T>(fct);
      }

      // Calculates the extent and communicates it to rank 0.
      template <class Writer>
      void writeLocalPieceImpl (Writer const& writer) const
      {
        auto&& extent = this->extent();

    #if HAVE_MPI
        int sendFlag = 0;
        MPI_Status sendStatus;
        MPI_Test(&sendRequest_, &sendFlag, &sendStatus);

        if (sendFlag) {
          int rank = gridView_.comm().rank();
          if (rank != 0) {
            MPI_Isend(extent.data(), extent.size(), MPI_INT, 0, /*tag=*/6, gridView_.comm(), &sendRequest_);
          } else {
            extents_[0] = extent;
          }
        }
    #endif

        writer(extent);
      }

      // Receive extent from all ranks and call the `writer` with the rank's extent vector
      template <class Writer>
      void writePiecesImpl (Writer const& writer) const
      {
    #if HAVE_MPI
        writer(0, extents_[0], true);

        int numRanks = gridView_.comm().size();
        for (int p = 1; p < numRanks; ++p) {
          int idx = -1;
          MPI_Status status;
          MPI_Waitany(numRanks, requests_.data(), &idx, &status);
          if (idx != MPI_UNDEFINED) {
            assert(idx == status.MPI_SOURCE && status.MPI_TAG == 6);
            writer(idx, extents_[idx], true);
          }
        }
    #else
        writer(0, this->extent(), true);
    #endif
      }

      // Origin (0,0,0)
      FieldVector<ctype, 3> originImpl () const
      {
        FieldVector<ctype, 3> vec; vec = ctype(0);
        return vec;
      }

      // Grid spacing (0,0,0)
      FieldVector<ctype, 3> spacingImpl () const
      {
        FieldVector<ctype, 3> vec; vec = ctype(0);
        return vec;
      }

      // Ordinate along each axis with constant \ref spacing from the \ref origin
      template <class T>
      std::array<std::vector<T>, 3> coordinatesImpl () const
      {
        auto origin = this->origin();
        auto spacing = this->spacing();
        auto extent = this->extent();

        std::array<std::vector<T>, 3> ordinates{};
        for (int d = 0; d < dim; ++d) {
          auto s = extent[2*d+1] - extent[2*d] + 1;
          ordinates[d].resize(s);
          for (int i = 0; i < s; ++i)
            ordinates[d][i] = origin[d] + (extent[2*d] + i)*spacing[d];
        }

        for (int d = dim; d < 3; ++d)
          ordinates[d].resize(1, T(0));

        return ordinates;
      }

    protected:
      using Super::gridView_;
      SubDataCollector subDataCollector_;

    #if HAVE_MPI
      mutable std::vector<std::array<int,6>> extents_;
      mutable std::vector<MPI_Request> requests_;
      mutable MPI_Request sendRequest_ = MPI_REQUEST_NULL;
    #endif
    };

    namespace Impl
    {
      // Should be specialized for concrete structured grid
      template <class GridView, class Grid>
      struct StructuredDataCollectorImpl;
    }

    template <class GridView>
    using StructuredDataCollector
      = typename Impl::StructuredDataCollectorImpl<GridView, typename GridView::Grid>::type;

  } // end namespace Vtk
} // end namespace Dune
