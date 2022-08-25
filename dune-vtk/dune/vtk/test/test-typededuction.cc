// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#if HAVE_CONFIG_H
#include "config.h" // autoconf defines, needed by the dune headers
#endif

#include <dune/vtk/vtkreader.hh>
#include <dune/vtk/vtkwriter.hh>

#if HAVE_DUNE_UGGRID
  #include <dune/grid/uggrid.hh>
  using GridType = Dune::UGGrid<2>;
#else
  #include <dune/grid/yaspgrid.hh>
  using GridType = Dune::YaspGrid<2>;
#endif
#include <dune/grid/utility/structuredgridfactory.hh>

int main (int argc, char** argv)
{
  using namespace Dune;
  MPIHelper::instance(argc, argv);

  auto grid = StructuredGridFactory<GridType>::createCubeGrid({0.0,0.0}, {1.0,2.0}, {2u,4u});

  // 1. construct writer from gridView
  VtkUnstructuredGridWriter writer1(grid->leafGridView());

  // 2. construct writer from datacollector
  Vtk::ContinuousDataCollector dataCollector1(grid->leafGridView());
  VtkUnstructuredGridWriter writer2(dataCollector1);
  VtkUnstructuredGridWriter writer3(stackobject_to_shared_ptr(dataCollector1));

  // 3. construct a default VtkWriter
  VtkWriter writer4(grid->leafGridView());

  // 4. construct reader from grid-factory
  GridFactory<GridType> factory;
  VtkReader reader1(factory);

  // 5. construct reader from grid-creator
  Vtk::ContinuousGridCreator creator(factory);
  VtkReader reader2(creator);
}