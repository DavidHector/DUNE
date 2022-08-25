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

#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/vtk/vtkreader.hh>
#include <dune/vtk/gridcreators/continuousgridcreator.hh>
#include <dune/vtk/gridcreators/discontinuousgridcreator.hh>
#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>

using namespace Dune;

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  const int dim = 3;

  using GridType = UGGrid<dim>;
  using GridView = typename GridType::LeafGridView;
  {
    FieldVector<double,dim> lowerLeft; lowerLeft = 0.0;
    FieldVector<double,dim> upperRight; upperRight = 1.0;
    auto numElements = filledArray<dim,unsigned int>(4);
    auto gridPtr = StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft, upperRight, numElements);
    auto& grid = *gridPtr;

    GridView gridView = grid.leafGridView();

    VtkUnstructuredGridWriter<GridView> vtkWriter1(gridView, Vtk::FormatTypes::ASCII);
    vtkWriter1.write("test_ascii_float32.vtu");

    VtkUnstructuredGridWriter<GridView> vtkWriter2(gridView, Vtk::FormatTypes::BINARY);
    vtkWriter2.write("test_binary_float32.vtu");

    VtkUnstructuredGridWriter<GridView> vtkWriter3(gridView, Vtk::FormatTypes::COMPRESSED, Vtk::DataTypes::FLOAT64);
    vtkWriter3.write("test_compressed_float64.vtu");
  }

  {
    std::cout << "read 'test_ascii_float32.vtu'..." << std::endl;
    auto gridPtr = VtkReader<GridType>::createGridFromFile("test_ascii_float32.vtu");
    auto& grid = *gridPtr;

    VtkUnstructuredGridWriter<GridView> vtkWriter(grid.leafGridView(), Vtk::FormatTypes::ASCII);
    vtkWriter.write("test_ascii_float32_2.vtu");
  }

  {
    std::cout << "read 'test_binary_float32.vtu'..." << std::endl;
    auto gridPtr = VtkReader<GridType>::createGridFromFile("test_binary_float32.vtu");
    auto& grid = *gridPtr;

    VtkUnstructuredGridWriter<GridView> vtkWriter(grid.leafGridView(), Vtk::FormatTypes::ASCII);
    vtkWriter.write("test_ascii_float32_3.vtu");
  }

  {
    std::cout << "read 'test_compressed_float64.vtu'..." << std::endl;
    auto gridPtr = VtkReader<GridType>::createGridFromFile("test_compressed_float64.vtu");
    auto& grid = *gridPtr;

    VtkUnstructuredGridWriter<GridView> vtkWriter(grid.leafGridView(), Vtk::FormatTypes::ASCII, Vtk::DataTypes::FLOAT64);
    vtkWriter.write("test_ascii_float64_3.vtu");
  }

  {
    std::cout << "read 'test_ascii_float32.vtu'..." << std::endl;
    auto gridPtr = VtkReader<GridType,Vtk::DiscontinuousGridCreator<GridType>>::createGridFromFile("test_ascii_float32.vtu");
    auto& grid = *gridPtr;

    VtkUnstructuredGridWriter<GridView> vtkWriter(grid.leafGridView(), Vtk::FormatTypes::ASCII);
    vtkWriter.write("test_ascii_float32_4.vtu");
  }

  if (Vtk::exists("paraview_3d.vtu")) {
    using GridType3d = UGGrid<3>;
    using GridView3d = typename GridType3d::LeafGridView;

    std::cout << "read 'paraview_3d.vtu'..." << std::endl;
    auto gridPtr = VtkReader<GridType3d>::createGridFromFile("paraview_3d.vtu");
    auto& grid = *gridPtr;

    VtkUnstructuredGridWriter<GridView3d> vtkWriter(grid.leafGridView(), Vtk::FormatTypes::ASCII, Vtk::DataTypes::FLOAT64);
    vtkWriter.write("paraview_3d_ascii.vtu");
  }

  if (Vtk::exists("paraview_2d.vtu")) {
    std::cout << "paraview_2d_ascii...\n";
    using GridType2d = UGGrid<2>;
    using GridView2d = typename GridType2d::LeafGridView;
    std::cout << "read 'paraview_2d.vtu'..." << std::endl;
    auto gridPtr = VtkReader<GridType2d>::createGridFromFile("paraview_2d.vtu");
    auto& grid = *gridPtr;

    VtkUnstructuredGridWriter<GridView2d> vtkWriter(grid.leafGridView(), Vtk::FormatTypes::ASCII, Vtk::DataTypes::FLOAT64);
    vtkWriter.write("paraview_2d_ascii.vtu");
  }
}
