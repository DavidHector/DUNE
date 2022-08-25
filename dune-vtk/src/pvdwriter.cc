// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <vector>

#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/filledarray.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/vtk/pvdwriter.hh>
#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>

using namespace Dune;
using namespace Dune::Functions;

template <class GridView>
void write (std::string prefix, GridView const& gridView)
{
  FieldVector<double,GridView::dimensionworld> c;
  if (GridView::dimensionworld > 0) c[0] = 11.0;
  if (GridView::dimensionworld > 1) c[1] = 7.0;
  if (GridView::dimensionworld > 2) c[2] = 3.0;

  // write analytic function
  auto p1Analytic = makeAnalyticGridViewFunction([&c](auto const& x) { return c.dot(x); }, gridView);

  using Writer = VtkUnstructuredGridWriter<GridView>;
  PvdWriter<Writer> pvdWriter(gridView, Vtk::FormatTypes::ASCII, Vtk::DataTypes::FLOAT32);

  std::string filename = prefix + "_" + std::to_string(GridView::dimensionworld) + "d_ascii.vtu";

  pvdWriter.addPointData(p1Analytic, "p1");
  pvdWriter.addCellData(p1Analytic, "p0");
  for (double t = 0.0; t < 10.0; t += 1.0) {
    pvdWriter.writeTimestep(t, filename);
  }
}

template <int I>
using int_ = std::integral_constant<int,I>;

int main (int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  // Test PvdWriter for YaspGrid
  using GridType = YaspGrid<2>;
  FieldVector<double,2> upperRight; upperRight = 1.0;
  auto numElements = filledArray<2,int>(8);
  GridType grid(upperRight, numElements, 0, 0);
  write("pvdwriter_yasp", grid.leafGridView());
}