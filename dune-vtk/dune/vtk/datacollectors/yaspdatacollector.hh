#pragma once

#include <array>
#include <type_traits>
#include <vector>

#include <dune/common/filledarray.hh>
#include <dune/common/fvector.hh>
#include <dune/common/std/type_traits.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/yaspgrid.hh>

#include "structureddatacollector.hh"

namespace Dune
{
  namespace Vtk
  {
    // Specialization for YaspGrid
    template <class GridView>
    class YaspDataCollector
        : public StructuredDataCollectorInterface<GridView, YaspDataCollector<GridView>>
    {
      using Self = YaspDataCollector;
      using Super = StructuredDataCollectorInterface<GridView, Self>;
      using ctype = typename GridView::ctype;

    public:
      using Super::dim;
      using Super::partition;

    public:
      YaspDataCollector (GridView const& gridView)
        : Super(gridView)
        , wholeExtent_(filledArray<6,int>(0))
        , extent_(filledArray<6,int>(0))
        , origin_(0.0)
        , spacing_(0.0)
        , level_(0)
      {}

      std::array<int, 6> const& wholeExtentImpl () const
      {
        return wholeExtent_;
      }

      std::array<int, 6> const& extentImpl () const
      {
        return extent_;
      }

      auto const& originImpl () const
      {
        return origin_;
      }

      auto const& spacingImpl () const
      {
        return spacing_;
      }

      void updateImpl ()
      {
        Super::updateImpl();

        level_ = gridView_.template begin<0,All_Partition>()->level();
        for (int i = 0; i < dim; ++i) {
          wholeExtent_[2*i] = 0;
          wholeExtent_[2*i+1] = grid(gridView_.grid()).levelSize(level_,i);
        }

        auto const& gl = *grid(gridView_.grid()).begin(level_);
        auto const& g = gl.interior[0];
        auto const& gc = *g.dataBegin();
        for (int i = 0; i < dim; ++i) {
          extent_[2*i] = gc.min(i);
          extent_[2*i+1] = gc.max(i)+1;
        }

        auto it = grid(gridView_.grid()).begin(level_);
        initGeometry(it->coords);
      }

      void initGeometry (EquidistantCoordinates<ctype,dim> const& coords)
      {
        for (int i = 0; i < dim; ++i) {
          spacing_[i] = coords.meshsize(i,0);
          origin_[i] = 0;
        }
      }

      void initGeometry (EquidistantOffsetCoordinates<ctype,dim> const& coords)
      {
        for (int i = 0; i < dim; ++i) {
          spacing_[i] = coords.meshsize(i,0);
          origin_[i] = coords.origin(i);
        }
      }

      void initGeometry (TensorProductCoordinates<ctype,dim> const& coords)
      {
        for (int i = 0; i < dim; ++i) {
          spacing_[i] = coords.meshsize(i,0); // is not constant, but also not used.
          origin_[i] = coords.coordinate(i,0);
        }
      }


      /// Extract the ordinates from the coordinates object of the current level
      template <class T>
      std::array<std::vector<T>, 3> coordinatesImpl () const
      {
        auto it = grid(gridView_.grid()).begin(level_);
        auto const& coords = it->coords;

        std::array<std::vector<T>, 3> ordinates{};
        for (int d = 0; d < dim; ++d) {
          auto s = extent_[2*d+1] - extent_[2*d] + 1;
          ordinates[d].resize(s);
          for (int i = 0; i < s; ++i)
            ordinates[d][i] = coords.coordinate(d, extent_[2*d] + i);
        }

        for (int d = dim; d < 3; ++d)
          ordinates[d].resize(1, T(0));

        return ordinates;
      }


    private:

      template <class G>
      using HostGrid = decltype(std::declval<G>().hostGrid());

      template <class G,
        std::enable_if_t<not Std::is_detected<HostGrid, G>::value, int> = 0>
      auto const& grid (G const& g) const
      {
        return g;
      }

      template <class G,
        std::enable_if_t<Std::is_detected<HostGrid, G>::value, int> = 0>
      auto const& grid (G const& g) const
      {
        return grid(g.hostGrid());
      }

    protected:
      using Super::gridView_;
      std::array<int, 6> wholeExtent_;
      std::array<int, 6> extent_;
      FieldVector<ctype,3> origin_;
      FieldVector<ctype,3> spacing_;
      int level_;
    };

    namespace Impl
    {
      template <class GridView, int dim, class Coordinates>
      struct StructuredDataCollectorImpl<GridView, YaspGrid<dim,Coordinates>>
      {
        using type = YaspDataCollector<GridView>;
      };
    }

  } // end namespace Vtk
} // end namespace Dune
