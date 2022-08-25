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
#include <dune/common/test/testsuite.hh>

// #include <dune/grid/uggrid.hh>
// #include <dune/alugrid/grid.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/vtk/vtkreader.hh>
#include <dune/vtk/datacollectors/lagrangedatacollector.hh>
#include <dune/vtk/gridcreators/lagrangegridcreator.hh>
#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>

using namespace Dune;

#ifndef GRID_PATH
#define GRID_PATH
#endif

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  const int dim = 2;
  const int dow = 3;
  const int order = 3;

  // using GridType = UGGrid<dim>;
  // using GridType = Dune::ALUGrid<dim,dow,Dune::simplex,Dune::conforming>;
  using GridType = FoamGrid<dim,dow>;
  using GridView = typename GridType::LeafGridView;
  using DataCollector = Vtk::LagrangeDataCollector<GridView, order>;
  using GridCreator = Vtk::LagrangeGridCreator<GridType>;

  std::string filename = "triangles_" + std::to_string(dow) + "d_order" + std::to_string(order);

  TestSuite testSuite;

  { // Test using the (static) file-reader interface
    std::cout << "Test 1..." << std::endl;
    auto gridPtr = VtkReader<GridType,GridCreator>::createGridFromFile(GRID_PATH "/" + filename + ".vtu");
    auto& grid = *gridPtr;
    grid.globalRefine(2);

    VtkUnstructuredGridWriter<GridView,DataCollector> vtkWriter(grid.leafGridView(), Vtk::FormatTypes::ASCII);
    vtkWriter.write(filename + "_out1.vtu");
  }

  { // Test using an instantiated reader
    std::cout << "Test 2..." << std::endl;
    GridFactory<GridType> factory;
    VtkReader<GridType,GridCreator> reader(factory);
    reader.read(GRID_PATH "/" + filename + ".vtu");
    auto gridPtr = factory.createGrid();
    auto& grid = *gridPtr;

    auto&& creator = reader.gridCreator();
    testSuite.check(creator.order() == order, "order");

    VtkUnstructuredGridWriter<GridView,DataCollector> vtkWriter(grid.leafGridView(), Vtk::FormatTypes::ASCII);
    vtkWriter.addPointData(creator, "param", 3);
    vtkWriter.write(filename + "_out2.vtu");
  }

  { // Test using an explicit grid-creator
    std::cout << "Test 3..." << std::endl;
    GridFactory<GridType> factory;
    GridCreator creator(factory);
    VtkReader<GridType,GridCreator> reader(creator);
    reader.read(GRID_PATH "/" + filename + ".vtu");
    auto gridPtr = factory.createGrid();
    auto& grid = *gridPtr;

    testSuite.check(creator.order() == order, "order");

    VtkUnstructuredGridWriter<GridView,DataCollector> vtkWriter(grid.leafGridView(), Vtk::FormatTypes::ASCII);
    vtkWriter.addPointData(creator, "param", 3);
    vtkWriter.write(filename + "_out3.vtu");
  }

  return testSuite.exit();
}
