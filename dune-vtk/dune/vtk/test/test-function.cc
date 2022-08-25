#include <config.h>

#include <optional>

#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/vtk/function.hh>
#include <dune/vtk/vtkwriter.hh>

#if HAVE_DUNE_UGGRID
  #include <dune/grid/uggrid.hh>
  using GridType = Dune::UGGrid<2>;
#else
  #include <dune/grid/yaspgrid.hh>
  using GridType = Dune::YaspGrid<2>;
#endif


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

int main (int argc, char** argv)
{
  using namespace Dune;
  MPIHelper::instance(argc, argv);

  auto grid = StructuredGridFactory<GridType>::createCubeGrid({0.0,0.0}, {1.0,2.0}, {2u,4u});

  using GridView = typename GridType::LeafGridView;
  GridView gridView = grid->leafGridView();
  VtkWriter writer{gridView, Vtk::FormatTypes::ASCII};

  // 1. add (legacy) VTKFunction
  std::vector<double> p1function(gridView.size(2), 1.0);
  using P1Function = P1VTKFunction<GridView,std::vector<double>>;
  std::shared_ptr<VTKFunction<GridView> const> f1(new P1Function(gridView, p1function, "p1"));
  writer.addPointData(f1);

  // 2. Add GlobalFunction
  auto f2 = GlobalFunction{gridView, [](auto x) { return x[0]+x[1]; }};
  writer.addPointData(f2, "global");

  // 3. Add data + Vtk::FieldInfo description
  auto f3 = GlobalFunction{gridView, [](auto x) { return x; }};
  writer.addPointData(f3, Dune::Vtk::FieldInfo{"vector1", 3, Vtk::RangeTypes::VECTOR});

  // 4. Add data + Dune::VTK::FieldInfo description
  writer.addPointData(f3, Dune::VTK::FieldInfo{"vector2", Dune::VTK::FieldInfo::Type::vector, 3, Dune::VTK::Precision::float64});

  // 5. Wrap function explicitly into Vtk::Function
  writer.addPointData(Vtk::Function<GridView>{f3, "vector3", 3}); // calls copy-constructor of Vtk::Function
  writer.addPointData(Vtk::Function<GridView>{f3, "vector3", 3}, "vector4", 3); // calls grid-function constructor

  writer.addPointData(Vtk::Function<GridView>{f3, "vector5a", {0,1}});
  writer.addPointData(Vtk::Function<GridView>{f3, "vector5"}, "vector5b", std::vector{0}, Vtk::RangeTypes::VECTOR);

  // 6. pass argument to FieldInfo and Function in any order
  writer.addPointData(f3, Dune::Vtk::FieldInfo{"vector6", Vtk::DataTypes::FLOAT32, Vtk::RangeTypes::VECTOR, 3});
  writer.addPointData(f3, "vector7", Vtk::DataTypes::FLOAT32, Vtk::RangeTypes::UNSPECIFIED, 3u);

  writer.addPointData(f3, "vector8", std::vector{0,1});

  // test default constructible
  Vtk::Function<GridView> func0;

  // test constructible by (legacy) VTKFunction
  Vtk::Function<GridView> func1{f1};
  VTK_ASSERT(func1.numComponents() == 1);
  VTK_ASSERT(func1.rangeType() == Vtk::RangeTypes::SCALAR);
  writer.addPointData(func1, "p2");

  // test constructible by local-function
  Vtk::Function<GridView> func2{f2, "func2"};
  VTK_ASSERT(func2.numComponents() == 1);
  VTK_ASSERT(func2.rangeType() == Vtk::RangeTypes::UNSPECIFIED);
  VTK_ASSERT(func2.dataType() == Vtk::DataTypes::FLOAT64);
  writer.addPointData(func2);

  Vtk::Function<GridView> func3{f3, "func3", 3};
  VTK_ASSERT(func3.numComponents() == 3);
  VTK_ASSERT(func3.rangeType() == Vtk::RangeTypes::UNSPECIFIED);
  VTK_ASSERT(func3.dataType() == Vtk::DataTypes::FLOAT64);
  writer.addPointData(func3);

  // test constructible with component vector
  Vtk::Function<GridView> func4a{f3, "func4a", 1};
  VTK_ASSERT(func4a.numComponents() == 1);
  VTK_ASSERT(func4a.rangeType() == Vtk::RangeTypes::UNSPECIFIED);
  VTK_ASSERT(func4a.dataType() == Vtk::DataTypes::FLOAT64);
  writer.addPointData(func4a);

  Vtk::Function<GridView> func4b{f3, "func4b", {1}}; // == func4a
  VTK_ASSERT(func4b.numComponents() == 1);
  VTK_ASSERT(func4b.rangeType() == Vtk::RangeTypes::UNSPECIFIED);
  VTK_ASSERT(func4b.dataType() == Vtk::DataTypes::FLOAT64);
  writer.addPointData(func4b);

  Vtk::Function<GridView> func4c{f3, "func4c", std::vector{1}};
  VTK_ASSERT(func4c.numComponents() == 1);
  VTK_ASSERT(func4c.rangeType() == Vtk::RangeTypes::UNSPECIFIED);
  VTK_ASSERT(func4c.dataType() == Vtk::DataTypes::FLOAT64);
  writer.addPointData(func4c);

  // Test copy-constructible
  auto func5{func3};
  VTK_ASSERT(func5.numComponents() == func3.numComponents());
  VTK_ASSERT(func5.rangeType() == func3.rangeType());
  VTK_ASSERT(func5.dataType() == func3.dataType());

  auto func6 = func3;
  VTK_ASSERT(func6.numComponents() == func3.numComponents());
  VTK_ASSERT(func6.rangeType() == func3.rangeType());
  VTK_ASSERT(func6.dataType() == func3.dataType());

  // Test move-constructible
  auto func7{std::move(func5)};
  VTK_ASSERT(func7.numComponents() == func3.numComponents());
  VTK_ASSERT(func7.rangeType() == func3.rangeType());
  VTK_ASSERT(func7.dataType() == func3.dataType());
  writer.addPointData(func7, "func7");

  auto func8 = std::move(func6);
  VTK_ASSERT(func8.numComponents() == func3.numComponents());
  VTK_ASSERT(func8.rangeType() == func3.rangeType());
  VTK_ASSERT(func8.dataType() == func3.dataType());
  writer.addPointData(func8, "func8");

  // test template-argument deduction
  Vtk::Function func9a{func8};
  VTK_ASSERT(func9a.numComponents() == func8.numComponents());
  VTK_ASSERT(func9a.rangeType() == func8.rangeType());
  VTK_ASSERT(func9a.dataType() == func8.dataType());

  Vtk::Function func9b{func8, "func9"};
  VTK_ASSERT(func9b.numComponents() == func8.numComponents());
  VTK_ASSERT(func9b.rangeType() == func8.rangeType());
  VTK_ASSERT(func9b.dataType() == func8.dataType());

  Vtk::Function func9c{func8, "func9", 2};
  VTK_ASSERT(func9c.numComponents() == 2);
  VTK_ASSERT(func9c.rangeType() == func8.rangeType());
  VTK_ASSERT(func9c.dataType() == func8.dataType());

  Vtk::Function func9d{func8, "func9", std::vector{0}};
  VTK_ASSERT(func9d.numComponents() == 1);
  VTK_ASSERT(func9d.rangeType() == func8.rangeType());
  VTK_ASSERT(func9d.dataType() == func8.dataType());

  Vtk::Function func9e{func8, "func9", Vtk::RangeTypes::SCALAR};
  VTK_ASSERT(func9e.numComponents() == 1);
  VTK_ASSERT(func9e.rangeType() == Vtk::RangeTypes::SCALAR);
  VTK_ASSERT(func9e.dataType() == func8.dataType());

  writer.write("test-function.vtu");
}