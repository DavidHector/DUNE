#pragma once

#include <type_traits>
#include <dune/geometry/type.hh>

namespace Dune
{
  namespace Vtk
  {
    template <class...> struct CheckTypes {};

    template <class DataCollector, class DC = std::decay_t<DataCollector>>
    using IsDataCollector = decltype((
      std::declval<DC&>().update(),
      std::declval<DC>().numPoints(),
      std::declval<DC>().numCells(),
      CheckTypes<typename DC::GridView>{},
    true));

    template <class GridView, class GV = std::decay_t<GridView>>
    using IsGridView = decltype((
      std::declval<GV>().grid(),
      std::declval<GV>().indexSet(),
      std::declval<GV>().size(0),
      std::declval<GV>().size(std::declval<Dune::GeometryType>()),
      CheckTypes<typename GV::Grid, typename GV::IndexSet>{},
    true));

    template <class GridFunction, class GF = std::decay_t<GridFunction>>
    using IsGridFunction = decltype((
      localFunction(std::declval<GF const&>()),
    true));

    template <class LocalFunction, class LocalContext, class LF = std::decay_t<LocalFunction>>
    using IsLocalFunction = decltype((
      std::declval<LF&>().bind(std::declval<LocalContext>()),
      std::declval<LF&>().unbind(),
      std::declval<LF>()(std::declval<typename LocalContext::Geometry::LocalCoordinate>()),
    true));

  } // end namespace Vtk
} // end namespace Dune
