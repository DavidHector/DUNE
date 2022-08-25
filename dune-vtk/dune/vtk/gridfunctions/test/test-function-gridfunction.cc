#include <config.h>

#if HAVE_DUNE_UGGRID
  #include <dune/grid/uggrid.hh>
  using GridType = Dune::UGGrid<2>;
#else
  #include <dune/grid/yaspgrid.hh>
  using GridType = Dune::YaspGrid<2>;
#endif

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/vtk/vtkreader.hh>
#include <dune/vtk/vtkwriter.hh>
#include <dune/vtk/function.hh>
#include <dune/vtk/gridcreators/lagrangegridcreator.hh>
#include <dune/vtk/utility/errors.hh>


// Wrapper for global-coordinate functions F
template <class GridView, class F>
class GlobalFunction
{
  using Element = typename GridView::template Codim<0>::Entity;
  using Geometry = typename Element::Geometry;

public:
  GlobalFunction (GridView const& gridView, F const& f)
    : gridView_(gridView)
    , f_(f)
  {}

  void bind(Element const& element) { geometry_.emplace(element.geometry()); }
  void unbind() { geometry_.reset(); }

  auto operator() (typename Geometry::LocalCoordinate const& local) const
  {
    assert(!!geometry_);
    return f_(geometry_->global(local));
  }

private:
  GridView gridView_;
  F f_;
  std::optional<Geometry> geometry_;
};

int main(int argc, char** argv)
{
  using namespace Dune;
  MPIHelper::instance(argc, argv);

  { // write point and cell data to a file

    auto grid = StructuredGridFactory<GridType>::createCubeGrid({0.0,0.0}, {1.0,2.0}, {2u,4u});
    auto gridView = grid->leafGridView();

    VtkUnstructuredGridWriter writer(gridView);

    auto f = GlobalFunction{gridView, [](auto x) { return x[0] + x[1]; }};
    writer.addPointData(f, "pointdata");
    writer.addCellData(f, "celldata");
    writer.write("test-function-gridfunction.vtu");
  }

  { // read data from file

    Vtk::LagrangeGridCreator<GridType> gridCreator;
    VtkReader reader(gridCreator);
    reader.read("test-function-gridfunction.vtu");

    auto grid = reader.createGrid();
    auto gridView = grid->leafGridView();
    using GridView = decltype(gridView);

    VtkUnstructuredGridWriter writer(gridView);

    // store point-data and cell-data gridfunction into a Vtk::Function
    Vtk::Function<GridView> pd{ reader.getPointData("pointdata") };
    Vtk::Function<GridView> cd{ reader.getCellData("celldata") };

    // store the grid-parametrization grid-function into a Vtk::Function
    Vtk::Function<GridView> param{ gridCreator };

    writer.addPointData(reader.getPointData("pointdata"));
    writer.addCellData(reader.getCellData("celldata"));
    writer.addPointData(gridCreator);

    auto localParam = localFunction(param);
    auto localPd = localFunction(pd);
    auto localCd = localFunction(cd);
    for (auto const& element : elements(gridView)) {
      localParam.bind(element);
      localPd.bind(element);
      localCd.bind(element);

      auto refElem = Dune::referenceElement(element);

      // evaluate grid-functions in local coordinates
      auto x0 = localParam(refElem.position(0,0));
      VTK_ASSERT(x0.size() == 2);

      auto x1 = localPd(refElem.position(0,0));
      VTK_ASSERT(x1.size() == 1);

      auto x2 = localCd(refElem.position(0,0));
      VTK_ASSERT(x2.size() == 1);
    }
  }
}