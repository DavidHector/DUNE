// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/common/hybridutilities.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_UGGRID
  #include <dune/grid/uggrid.hh>
#endif

#include <dune/vtk/datacollectors/continuousdatacollector.hh>
#include <dune/vtk/datacollectors/lagrangedatacollector.hh>

/**
 * This test compares LagrangeDataCollector<1> against the ContinuousDataCollector, thus
 * just compares the  linear order implementations
 **/

template <class DataCollector1, class DataCollector2>
void testDataCollector (std::string prefix, Dune::TestSuite& testSuite, DataCollector1& dataCollector1, DataCollector2& dataCollector2)
{
  using std::sqrt;
  auto tol = sqrt(std::numeric_limits<double>::epsilon());

  dataCollector1.update();
  dataCollector2.update();

  // check that number of points and cells are equal
  testSuite.check(dataCollector1.numPoints() == dataCollector2.numPoints(), prefix + "_numPoints");
  testSuite.check(dataCollector1.numCells() == dataCollector2.numCells(), prefix + "_numCells");

  auto points1 = dataCollector1.template points<double>();
  auto points2 = dataCollector2.template points<double>();

  // check that point sizes are equal
  testSuite.check(points1.size() == points2.size(), prefix + "_points.size");
  testSuite.check(points1.size() == 3*dataCollector1.numPoints(), prefix + "_points.size/3");

  // check that point coordinates are equal
  using std::abs;
  for (std::size_t i = 0; i < points1.size(); ++i)
    testSuite.check(abs(points1[i] - points2[i]) < tol, prefix + "_points[" + std::to_string(i) + "]");

  auto cells1 = dataCollector1.cells();
  auto cells2 = dataCollector2.cells();

  // check that cell sizes are equal
  testSuite.check(cells1.types.size() == cells2.types.size(), prefix + "_cells.types.size");
  testSuite.check(cells1.offsets.size() == cells2.offsets.size(), prefix + "_cells.offsets.size");
  testSuite.check(cells1.connectivity.size() == cells2.connectivity.size(), prefix + "_cells.connectivity.size");

  // NOTE: cells.types do not need to be equal, e.g. LINEAR != LAGRANGE_LINEAR

  // check that offsets are equal
  for (std::size_t i = 0; i < cells1.offsets.size(); ++i)
    testSuite.check(cells1.offsets[i] == cells2.offsets[i], prefix + "_cells.offsets[" + std::to_string(i) + "]");

  // check that connectivities are equal
  for (std::size_t i = 0; i < cells1.connectivity.size(); ++i)
    testSuite.check(cells1.connectivity[i] == cells2.connectivity[i], prefix + "_cells.connectivity[" + std::to_string(i) + "]");
}

template <class GridView>
void testGridView (std::string prefix, Dune::TestSuite& testSuite, GridView const& gridView)
{
  // 1. test linear order lagrange data-collector
  {
    Dune::Vtk::ContinuousDataCollector<GridView> linearDataCollector(gridView);

    // data collector with template order
    Dune::Vtk::LagrangeDataCollector<GridView, 1> lagrangeDataCollector1a(gridView);
    testDataCollector(prefix + "_linear_template", testSuite, lagrangeDataCollector1a, linearDataCollector);

    // data collector with runtime order
    Dune::Vtk::LagrangeDataCollector<GridView> lagrangeDataCollector1b(gridView, 1);
    testDataCollector(prefix + "_linear_runtime", testSuite, lagrangeDataCollector1b, linearDataCollector);
  }
}

template <class Grid>
void testGrid (std::string prefix, Dune::TestSuite& testSuite, Grid& grid)
{
  grid.globalRefine(1);

  testGridView(prefix + "_level0", testSuite, grid.levelGridView(0));
  testGridView(prefix + "_leaf", testSuite, grid.leafGridView());
}

int main(int argc, char** argv)
{
  using namespace Dune;
  MPIHelper::instance(argc, argv);

  TestSuite testSuite;

  YaspGrid<1> yaspGrid1({1.0}, {1});
  YaspGrid<2> yaspGrid2({1.0,1.0}, {1,1});
  YaspGrid<3> yaspGrid3({1.0,1.0,1.0}, {1,1,1});

  testGrid("yaspgrid_1d", testSuite, yaspGrid1);
  testGrid("yaspgrid_2d", testSuite, yaspGrid2);
  testGrid("yaspgrid_3d", testSuite, yaspGrid3);

#if HAVE_DUNE_UGGRID
  auto ugGrid2 = StructuredGridFactory<UGGrid<2>>::createSimplexGrid({0.0,0.0}, {1.0,1.0}, {1u,1u});
  auto ugGrid3 = StructuredGridFactory<UGGrid<3>>::createSimplexGrid({0.0,0.0,0.0}, {1.0,1.0,1.0}, {1u,1u,1u});

  testGrid("uggrid_2d", testSuite, *ugGrid2);
  testGrid("uggrid_3d", testSuite, *ugGrid3);
#endif

  return testSuite.exit();
}
