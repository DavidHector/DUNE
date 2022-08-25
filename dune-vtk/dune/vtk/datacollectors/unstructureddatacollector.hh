#pragma once

#include <cstdint>
#include <vector>

#include <dune/grid/common/partitionset.hh>
#include <dune/vtk/datacollectorinterface.hh>

namespace Dune
{
  namespace Vtk
  {
    struct Cells
    {
      std::vector<std::uint8_t> types;
      std::vector<std::int64_t> offsets;
      std::vector<std::int64_t> connectivity;
    };

    template <class GridView, class Derived, class Partition = Partitions::InteriorBorder>
    class UnstructuredDataCollectorInterface
        : public DataCollectorInterface<GridView, Derived, Partition>
    {
      using Super = DataCollectorInterface<GridView, Derived, Partition>;

    public:
      using Super::dim;
      using Super::partition;

    public:
      UnstructuredDataCollectorInterface (GridView const& gridView)
        : Super(gridView)
      {}

      /// \brief Return cell types, offsets, and connectivity. \see Cells
      Cells cells () const
      {
        return this->asDerived().cellsImpl();
      }

      std::vector<std::uint64_t> pointIds () const
      {
        return this->asDerived().pointIdsImpl();
      }

    protected:
      // default implementation
      std::vector<std::uint64_t> pointIdsImpl () const
      {
        return {};
      }

    protected:
      using Super::gridView_;
    };

  } // end namespace Vtk
} // end namespace Dune
