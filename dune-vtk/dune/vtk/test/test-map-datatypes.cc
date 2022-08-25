#include <config.h>

#include <type_traits>
#include <dune/vtk/types.hh>


int main()
{
  using namespace Dune;

  // unconstrained mapping of DataTypes to real type t
  Vtk::mapDataTypes(Vtk::DataTypes::UINT32, [](auto t)
  {
    using T = typename decltype(t)::type;
    VTK_ASSERT(Vtk::dataTypeOf<T>() == Vtk::DataTypes::UINT32);
  });

  // constrained mapping of DataTypes to real type t
  Vtk::mapDataTypes<std::is_integral>(Vtk::DataTypes::UINT32, [](auto t)
  {
    using T = typename decltype(t)::type;
    VTK_ASSERT(Vtk::dataTypeOf<T>() == Vtk::DataTypes::UINT32);
    static_assert(std::is_integral_v<T>);
  });

  Vtk::mapDataTypes<std::is_floating_point>(Vtk::DataTypes::FLOAT32, [](auto t)
  {
    using T = typename decltype(t)::type;
    VTK_ASSERT(Vtk::dataTypeOf<T>() == Vtk::DataTypes::FLOAT32);
    static_assert(std::is_floating_point_v<T>);
  });


  // if the DataType does not fulfill the constraint, the function is never invoked
  Vtk::mapDataTypes<std::is_floating_point>(Vtk::DataTypes::UINT32, [](auto t)
  {
    // This function should never be called
    VTK_ASSERT(false);
  });


  // unconstrained mapping with multiple types
  Vtk::mapDataTypes(Vtk::DataTypes::FLOAT32, Vtk::DataTypes::UINT32,
  [](auto f, auto h)
  {
    using F = typename decltype(f)::type;
    using H = typename decltype(h)::type;
    VTK_ASSERT(Vtk::dataTypeOf<F>() == Vtk::DataTypes::FLOAT32);
    VTK_ASSERT(Vtk::dataTypeOf<H>() == Vtk::DataTypes::UINT32);
  });

  // constrained mapping with multiple types
  Vtk::mapDataTypes<std::is_floating_point,std::is_integral>(Vtk::DataTypes::FLOAT32, Vtk::DataTypes::UINT32,
  [](auto f, auto h)
  {
    using F = typename decltype(f)::type;
    using H = typename decltype(h)::type;
    VTK_ASSERT(Vtk::dataTypeOf<F>() == Vtk::DataTypes::FLOAT32);
    VTK_ASSERT(Vtk::dataTypeOf<H>() == Vtk::DataTypes::UINT32);
    static_assert(std::is_floating_point_v<F>);
    static_assert(std::is_integral_v<H>);
  });
}