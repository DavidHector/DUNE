// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <set>
#include <vector>

#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions

#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/vtk/vtkwriter.hh>

using namespace Dune;
using namespace Dune::Functions;

using TestCases = std::set<std::tuple<std::string,Vtk::FormatTypes,Vtk::DataTypes>>;
static TestCases test_cases = {
  {"ascii32", Vtk::FormatTypes::ASCII, Vtk::DataTypes::FLOAT32},
  {"bin32", Vtk::FormatTypes::BINARY, Vtk::DataTypes::FLOAT32},
  {"zlib32", Vtk::FormatTypes::COMPRESSED, Vtk::DataTypes::FLOAT32},
  {"ascii64", Vtk::FormatTypes::ASCII, Vtk::DataTypes::FLOAT64},
  {"bin64", Vtk::FormatTypes::BINARY, Vtk::DataTypes::FLOAT64},
  {"zlib64", Vtk::FormatTypes::COMPRESSED, Vtk::DataTypes::FLOAT64},
};

template <class GridView>
void write (std::string prefix, GridView const& gridView)
{
  FieldVector<double,GridView::dimensionworld> c;
  if (GridView::dimensionworld > 0) c[0] = 11.0;
  if (GridView::dimensionworld > 1) c[1] = 7.0;
  if (GridView::dimensionworld > 2) c[2] = 3.0;

  auto vec = [&c](auto const& x) {
    FieldVector<double,GridView::dimensionworld> result; result = c.dot(x);
    return result;
  };
  auto mat = [&c](auto const& x) {
    FieldMatrix<double,GridView::dimensionworld,GridView::dimensionworld> result;
    for (int i = 0; i < GridView::dimensionworld; ++i)
      result[i][i] = c.dot(x);
    return result;
  };
  auto vec_fct = makeAnalyticGridViewFunction(vec, gridView);
  auto mat_fct = makeAnalyticGridViewFunction(mat, gridView);

  for (auto const& test_case : test_cases) {
    VtkWriter<GridView> vtkWriter(gridView, std::get<1>(test_case), std::get<2>(test_case));
    vtkWriter.addPointData(vec_fct, "vec");
    vtkWriter.addPointData(mat_fct, "mat");
    vtkWriter.write(prefix + "_" + std::to_string(GridView::dimensionworld) + "d_" + std::get<0>(test_case) + ".vtu");
  }
}

template <int I>
using int_ = std::integral_constant<int,I>;

int main (int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  // Test VtkWriter for YaspGrid
  Hybrid::forEach(std::make_tuple(int_<2>{}, int_<3>{}), [](auto dim)
  {
    using GridType = YaspGrid<dim.value>;
    FieldVector<double,dim.value> upperRight; upperRight = 1.0;
    auto numElements = filledArray<dim.value,int>(8);
    GridType grid(upperRight, numElements, 0, 0);
    write("vectorwriter_yasp", grid.leafGridView());
  });
}