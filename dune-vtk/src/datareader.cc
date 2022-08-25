// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <cmath>
#include <iostream>
#include <vector>

#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/filledarray.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/vtk/vtkreader.hh>
#include <dune/vtk/datacollectors/lagrangedatacollector.hh>
#include <dune/vtk/gridcreators/continuousgridcreator.hh>
#include <dune/vtk/gridcreators/lagrangegridcreator.hh>
#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>

#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

using namespace Dune;

template <class GridType, class GridCreator>
void run(std::string const& prefix)
{
  std::tuple f_tuple{0,
    [](auto const& x) { return std::sin(x[0]*x[0]*x[0] + 20); },
    [](auto const& x) { return std::sin(x[0]*x[0]) + std::cos(x[1] + x[0]); },
    [](auto const& x) { return std::sin(x[0] + x[1] + x[2]) + std::cos(x[0]*x[1]*x[2]) + std::sin(x[0]*x[2])*std::cos(x[1] + x[2]); }
  };

  constexpr int dim = GridType::dimension;
  auto f = std::get<dim>(f_tuple);

  using GridView = typename GridType::LeafGridView;
  {
    FieldVector<double,dim> lowerLeft; lowerLeft = 0.0;
    FieldVector<double,dim> upperRight; upperRight = 1.0;
    auto numElements = filledArray<dim,unsigned int>(4);
    auto gridPtr = StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft, upperRight, numElements);
    auto& grid = *gridPtr;

    GridView gridView = grid.leafGridView();
    auto data = Functions::makeAnalyticGridViewFunction(f, gridView);

    VtkUnstructuredGridWriter<GridView> vtkWriter1(gridView, Vtk::FormatTypes::ASCII);
    vtkWriter1.addPointData(data, "f");
    vtkWriter1.addCellData(data, "g");
    vtkWriter1.write(prefix + "_ascii.vtu");

    VtkUnstructuredGridWriter<GridView> vtkWriter2(gridView, Vtk::FormatTypes::BINARY);
    vtkWriter2.addPointData(data, "f");
    vtkWriter2.addCellData(data, "g");
    vtkWriter2.write(prefix + "_binary.vtu");

    VtkUnstructuredGridWriter<GridView> vtkWriter3(gridView, Vtk::FormatTypes::COMPRESSED, Vtk::DataTypes::FLOAT64);
    vtkWriter3.addPointData(data, "f");
    vtkWriter3.addCellData(data, "g");
    vtkWriter3.write(prefix + "_compressed.vtu");

    if constexpr(std::is_same_v<GridCreator, Vtk::LagrangeGridCreator<GridType>>) {
      VtkUnstructuredGridWriter<GridView, Vtk::LagrangeDataCollector<GridView,3>> vtkWriter4(gridView, Vtk::FormatTypes::BINARY, Vtk::DataTypes::FLOAT64);
      vtkWriter4.addPointData(data, "f");
      vtkWriter4.addCellData(data, "g");
      vtkWriter4.write(prefix + "_order3.vtu");
    }
  }

  {
    std::cout << "read '" << prefix << "_ascii.vtu'..." << std::endl;
    VtkReader<GridType, GridCreator, float> reader;
    reader.read(prefix + "_ascii.vtu");
    auto gridPtr = reader.createGrid();
    auto& grid = *gridPtr;

    auto gf = reader.getPointData("f");
    auto gg = reader.getCellData("g");

    VtkUnstructuredGridWriter<GridView> vtkWriter(grid.leafGridView(), Vtk::FormatTypes::ASCII);
    vtkWriter.addPointData(gf, "f");
    vtkWriter.addCellData(gg, "g");
    vtkWriter.write(prefix + "_ascii_out.vtu");
  }

  {
    std::cout << "read '" << prefix << "_binary.vtu'..." << std::endl;
    VtkReader<GridType, GridCreator, float> reader;
    reader.read(prefix + "_binary.vtu");
    auto gridPtr = reader.createGrid();
    auto& grid = *gridPtr;

    auto gf = reader.getPointData("f");
    auto gg = reader.getCellData("g");

    VtkUnstructuredGridWriter<GridView> vtkWriter(grid.leafGridView(), Vtk::FormatTypes::ASCII);
    vtkWriter.addPointData(gf, "f");
    vtkWriter.addCellData(gg, "g");
    vtkWriter.write(prefix + "_binary_out.vtu");
  }

  {
    std::cout << "read '" << prefix << "_compressed.vtu'..." << std::endl;
    VtkReader<GridType, GridCreator> reader;
    reader.read(prefix + "_compressed.vtu");
    auto gridPtr = reader.createGrid();
    auto& grid = *gridPtr;

    auto gf = reader.getPointData("f");
    auto gg = reader.getCellData("g");

    VtkUnstructuredGridWriter<GridView> vtkWriter(grid.leafGridView(), Vtk::FormatTypes::ASCII);
    vtkWriter.addPointData(gf, "f");
    vtkWriter.addCellData(gg, "g");
    vtkWriter.write(prefix + "_compressed_out.vtu");
  }

  if constexpr(std::is_same_v<GridCreator, Vtk::LagrangeGridCreator<GridType>>) {
    std::cout << "read '" << prefix << "_order3.vtu'..." << std::endl;
    VtkReader<GridType, GridCreator> reader;
    reader.read(prefix + "_order3.vtu");
    auto gridPtr = reader.createGrid();
    auto& grid = *gridPtr;

    auto gf = reader.getPointData("f");
    auto gg = reader.getCellData("g");

    VtkUnstructuredGridWriter<GridView, Vtk::LagrangeDataCollector<GridView,3>> vtkWriter(grid.leafGridView(), Vtk::FormatTypes::ASCII);
    vtkWriter.addPointData(gf, "f");
    vtkWriter.addCellData(gg, "g");
    vtkWriter.write(prefix + "_order3_out.vtu");
  }
}

template <class GridType>
void run(std::string const& prefix)
{
  run<GridType, Vtk::ContinuousGridCreator<GridType>>(prefix + "_continuous");
  run<GridType, Vtk::LagrangeGridCreator<GridType>>(prefix + "_lagrange");
}

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  run<UGGrid<2>>("uggrid_2d");
  run<UGGrid<3>>("uggrid_3d");

#if HAVE_DUNE_ALUGRID
  run<Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>>("alugrid_2d");
  run<Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>>("alugrid_3d");
#endif
}
