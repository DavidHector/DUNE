#pragma once

#include <dune/vtk/writers/vtkimagedatawriter.hh>
#include <dune/vtk/writers/vtkrectilineargridwriter.hh>
#include <dune/vtk/writers/vtkstructuredgridwriter.hh>
#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>

#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid.hh>
#include <dune/vtk/datacollectors/spdatacollector.hh>
#endif

#include <dune/grid/geometrygrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/vtk/datacollectors/yaspdatacollector.hh>

namespace Dune
{
  namespace Impl
  {
    // The default writer assumes an unstructured grid
    template <class GridView, class Grid>
    struct VtkWriterImpl
    {
      using type = VtkUnstructuredGridWriter<GridView>;
    };
  }

  /// \brief Default choice of VTK Writer for several grid types.
  /**
   * Choose a VTK writer depending on the grid type. Some specialization for standard dune-grid grids
   * are provided, like YaspGrid and GeometrGrid.
   *
   * Note: Uses the default data-collector. If you want to choose a special data-collector, use
   * the concrete writer Implementation instead. \see VtkWriterInterface
   **/
  template <class GridView>
  class VtkWriter
      : public Impl::VtkWriterImpl<GridView, typename GridView::Grid>::type
  {
    using Super = typename Impl::VtkWriterImpl<GridView, typename GridView::Grid>::type;

  public:
    using Super::Super;
  };

  // deduction guide
  template <class GridView,
    class = std::void_t<typename GridView::IndexSet>>
  VtkWriter (GridView const&, Vtk::FormatTypes = Vtk::FormatTypes::BINARY, Vtk::DataTypes = Vtk::DataTypes::FLOAT32)
    -> VtkWriter<GridView>;


  namespace Impl
  {
    // A structured grid with coordinates in x, y, and z direction with arbitrary spacing
    template <class GridView, int dim, class Coordinates>
    struct VtkWriterImpl<GridView, YaspGrid<dim,Coordinates>>
    {
      using type = VtkRectilinearGridWriter<GridView, Vtk::YaspDataCollector<GridView>>;
    };

#if HAVE_DUNE_SPGRID
    // A structured grid with coordinates in x, y, and z direction with arbitrary spacing
    template <class GridView, class ct, int dim, template <int> class Ref, class Comm>
    struct VtkWriterImpl<GridView, SPGrid<ct,dim,Ref,Comm>>
    {
      using type = VtkRectilinearGridWriter<GridView, Vtk::SPDataCollector<GridView>>;
    };
#endif

    // A structured grid with coordinates in x, y, and z direction with arbitrary spacing
    template <class GridView, int dim, class ct>
    struct VtkWriterImpl<GridView, YaspGrid<dim,TensorProductCoordinates<ct,dim>>>
    {
      using type = VtkRectilinearGridWriter<GridView, Vtk::YaspDataCollector<GridView>>;
    };

    // A transformed structured grid has structured connectivity but unstructured point
    // coordinates.
    template <class GridView, int dim, class Coordinates, class CoordFunction, class Allocator>
    struct VtkWriterImpl<GridView, GeometryGrid<YaspGrid<dim,Coordinates>, CoordFunction, Allocator>>
    {
      using type = VtkStructuredGridWriter<GridView, Vtk::YaspDataCollector<GridView>>;
    };

  } // end namespace Impl
} // end namespace Dune
