// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#if HAVE_CONFIG_H
#include "config.h" // autoconf defines, needed by the dune headers
#endif

#include <algorithm>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_UGGRID
  #include <dune/grid/uggrid.hh>
  #include <dune/grid/utility/structuredgridfactory.hh>
#endif
#include <dune/vtk/vtkwriter.hh>

#include "checkvtkfile.hh"

// accumulate exit status
void acc (int &accresult, int result)
{
  if (accresult == 0 || (accresult == 77 && result != 0))
    accresult = result;
}

struct Acc
{
  int operator() (int v1, int v2) const
  {
    acc(v1, v2);
    return v1;
  }
};

template <class GridView>
int testGridView (std::string prefix, Dune::VTKChecker& vtkChecker, GridView const& gridView)
{
  enum { dim = GridView :: dimension };

  Dune::VtkWriter<GridView> vtk(gridView);

  auto f1 = Dune::Functions::makeAnalyticGridViewFunction([](const auto& x) { return std::sin(x.two_norm()); },gridView);
  vtk.addCellData(f1, "scalar-valued lambda");

  auto f2 = Dune::Functions::makeAnalyticGridViewFunction([](const auto& x) { return x; },gridView);
  vtk.addPointData(f2, "vector-valued lambda");

  int result = 0;

  // ASCII files
  {
    vtk.setFormat(Dune::Vtk::FormatTypes::ASCII);
    std::string name = vtk.write(prefix + "_ascii");
    if (gridView.comm().rank() == 0) vtkChecker.push(name);
  }

  // BINARY files
  {
    vtk.setFormat(Dune::Vtk::FormatTypes::BINARY);
    std::string name = vtk.write(prefix + "_binary");
    if (gridView.comm().rank() == 0) vtkChecker.push(name);
  }

  // COMPRESSED files
  {
    vtk.setFormat(Dune::Vtk::FormatTypes::COMPRESSED);
    std::string name = vtk.write(prefix + "_compressed");
    if (gridView.comm().rank() == 0) vtkChecker.push(name);
  }

  return result;
}

template <class Grid>
int testGrid (std::string prefix, Dune::VTKChecker& vtkChecker, Grid& grid)
{
  if (grid.comm().rank() == 0)
    std::cout << "vtkCheck(" << prefix << ")" << std::endl;

  grid.globalRefine(1);

  int result = 0;

  acc(result, testGridView( prefix + "_leaf",   vtkChecker, grid.leafGridView() ));
  acc(result, testGridView( prefix + "_level0", vtkChecker, grid.levelGridView(0) ));
  acc(result, testGridView( prefix + "_level1", vtkChecker, grid.levelGridView(grid.maxLevel()) ));

  return result;
}

int main (int argc, char** argv)
{
  using namespace Dune;
  const MPIHelper& mpiHelper = MPIHelper::instance(argc, argv);

  if (mpiHelper.rank() == 0)
    std::cout << "vtktest: MPI_Comm_size == " << mpiHelper.size() << std::endl;

  int result = 0; // pass by default

  VTKChecker vtkChecker;

  YaspGrid<1> yaspGrid1({1.0}, {8});
  YaspGrid<2> yaspGrid2({1.0,2.0}, {8,4});
  YaspGrid<3> yaspGrid3({1.0,2.0,3.0}, {8,4,4});

  acc(result, testGrid("yaspgrid_1d", vtkChecker, yaspGrid1));
  acc(result, testGrid("yaspgrid_2d", vtkChecker, yaspGrid2));
  acc(result, testGrid("yaspgrid_3d", vtkChecker, yaspGrid3));

#if HAVE_DUNE_UGGRID
  auto ugGrid2a = StructuredGridFactory<UGGrid<2>>::createSimplexGrid({0.0,0.0}, {1.0,2.0}, {2u,4u});
  auto ugGrid3a = StructuredGridFactory<UGGrid<3>>::createSimplexGrid({0.0,0.0,0.0}, {1.0,2.0,3.0}, {2u,4u,6u});

  acc(result, testGrid("uggrid_simplex_2d", vtkChecker, *ugGrid2a));
  acc(result, testGrid("uggrid_simplex_3d", vtkChecker, *ugGrid3a));

  auto ugGrid2b = StructuredGridFactory<UGGrid<2>>::createCubeGrid({0.0,0.0}, {1.0,2.0}, {2u,4u});
  auto ugGrid3b = StructuredGridFactory<UGGrid<3>>::createCubeGrid({0.0,0.0,0.0}, {1.0,2.0,3.0}, {2u,4u,6u});

  acc(result, testGrid("uggrid_cube_2d", vtkChecker, *ugGrid2b));
  acc(result, testGrid("uggrid_cube_3d", vtkChecker, *ugGrid3b));
#endif

  acc(result, vtkChecker.check());

  mpiHelper.getCollectiveCommunication().allreduce<Acc>(&result, 1);

  return result;
}
