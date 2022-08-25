// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <cstring>
#include <iostream>
#include <vector>

#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/filledarray.hh>
#include <dune/common/std/type_traits.hh>
#include <dune/common/test/testsuite.hh>

#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dune/vtk/vtkreader.hh>
#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>
#include <dune/vtk/gridcreators/parallelgridcreator.hh>
#include <dune/vtk/gridcreators/serialgridcreator.hh>

using namespace Dune;

// see https://stackoverflow.com/questions/6163611/compare-two-files
bool compare_files (std::string const& fn1, std::string const& fn2)
{
  std::ifstream in1(fn1, std::ios::binary);
  std::ifstream in2(fn2, std::ios::binary);
  if (!in1 || !in2) {
    std::cout << "can not find file " << fn1 << " or file " << fn2 << "\n";
    return false;
  }

  std::ifstream::pos_type size1 = in1.seekg(0, std::ifstream::end).tellg();
  in1.seekg(0, std::ifstream::beg);

  std::ifstream::pos_type size2 = in2.seekg(0, std::ifstream::end).tellg();
  in2.seekg(0, std::ifstream::beg);

  if (size1 != size2)
    return false;

  static const std::size_t BLOCKSIZE = 4096;
  std::size_t remaining = size1;

  while (remaining) {
    char buffer1[BLOCKSIZE], buffer2[BLOCKSIZE];
    std::size_t size = std::min(BLOCKSIZE, remaining);

    in1.read(buffer1, size);
    in2.read(buffer2, size);

    if (0 != std::memcmp(buffer1, buffer2, size))
      return false;

    remaining -= size;
  }

  return true;
}


template <class GF>
using HasParallelGridFactoryImpl = decltype(std::declval<GF>().createGrid(true,true,std::string(""),true));

template <class G>
using HasParallelGridFactory = Std::is_detected<HasParallelGridFactoryImpl, GridFactory<G>>;


template <class Test>
void compare (Test& test, Vtk::Path const& dir, Vtk::Path const& name)
{
  test.check(compare_files(dir.string() + '/' + name.string() + ".vtu",
                           dir.string() + '/' + name.string() + "_2.vtu"));
}

template <class GridView>
void writer_test (GridView const& gridView, std::string base_name)
{
  VtkUnstructuredGridWriter<GridView> vtkWriter(gridView, Vtk::FormatTypes::ASCII, Vtk::DataTypes::FLOAT32);
  vtkWriter.write(base_name + ".vtu");
}

template <class Grid, class Creator>
void reader_writer_test(MPIHelper& mpi, TestSuite& test, std::string const& testName, bool doLoadBalance = true)
{
  std::cout << "== " << testName << "\n";
  std::string base_name = "parallel_rw_dim" + std::to_string(Grid::dimension);
  std::vector<std::string> pieces1, pieces2;

  std::string ext = ".vtu";
  if (mpi.size() > 1)
    ext = ".pvtu";

  // Step 1: create a new grid and write it to file1
  const int dim = Grid::dimension;
  {
    FieldVector<double,dim> lowerLeft; lowerLeft = 0.0;
    FieldVector<double,dim> upperRight; upperRight = 1.0;
    auto numElements = filledArray<dim,unsigned int>(4);
    auto gridPtr = StructuredGridFactory<Grid>::createSimplexGrid(lowerLeft, upperRight, numElements);
    gridPtr->loadBalance();
    writer_test(gridPtr->leafGridView(), base_name);
  }

  mpi.getCollectiveCommunication().barrier(); // need a barrier between write and read

  // Step 2: read the grid from file1 and write it back to file2
  { GridFactory<Grid> factory;
    VtkReader<Grid, Creator> reader{factory};
    reader.read(base_name + ext);

    std::unique_ptr<Grid> grid{ Hybrid::ifElse(HasParallelGridFactory<Grid>{},
      [&](auto id) { return id(factory).createGrid(std::true_type{}); },
      [&](auto id) { return id(factory).createGrid(); }) };
    if (doLoadBalance)
      grid->loadBalance();

    writer_test(grid->leafGridView(), base_name + "_1");
  }

  mpi.getCollectiveCommunication().barrier();

  // Step 3: read the (parallel) file1 to get the piece filenames
  { GridFactory<Grid> factory;
    VtkReader<Grid, Creator> reader{factory};
    reader.read(base_name + "_1" + ext);

    std::unique_ptr<Grid> grid{ Hybrid::ifElse(HasParallelGridFactory<Grid>{},
      [&](auto id) { return id(factory).createGrid(std::true_type{}); },
      [&](auto id) { return id(factory).createGrid(); }) };
    if (doLoadBalance)
      grid->loadBalance();
    pieces1 = reader.pieces();

    writer_test(grid->leafGridView(), base_name + "_2");
  }

  mpi.getCollectiveCommunication().barrier();

  // Step 4: read the (parallel) file2 to get the piece filenames
  { GridFactory<Grid> factory;
    VtkReader<Grid, Creator> reader{factory};
    reader.read(base_name + "_2" + ext, false);

    pieces2 = reader.pieces();
  }

  mpi.getCollectiveCommunication().barrier();

  // Step 4: compare the pieces
  if (mpi.rank() == 0) {
    test.check(pieces1.size() == pieces2.size(), "pieces1.size == pieces2.size");
    for (std::size_t i = 0; i < pieces1.size(); ++i)
      test.check(compare_files(pieces1[i], pieces2[i]), "compare(" + pieces1[i] + ", " + pieces2[i] + ")");
  }
}

#if HAVE_DUNE_ALUGRID
  template <int dim>
  using ALUGridType = Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming>;
#endif

template <int I>
using int_ = std::integral_constant<int,I>;

int main (int argc, char** argv)
{
  auto& mpi = Dune::MPIHelper::instance(argc, argv);

  TestSuite test{};

#if HAVE_DUNE_UGGRID
  reader_writer_test<UGGrid<2>, Vtk::SerialGridCreator<UGGrid<2>>>(mpi, test, "UGGrid<2>");
  reader_writer_test<UGGrid<3>, Vtk::SerialGridCreator<UGGrid<3>>>(mpi, test, "UGGrid<3>");
#endif

#if HAVE_DUNE_ALUGRID
  // Test VtkWriter for ALUGrid.
  reader_writer_test<ALUGridType<2>, Vtk::SerialGridCreator<ALUGridType<2>>>(mpi, test, "ALUGridType<2>");
  reader_writer_test<ALUGridType<2>, Vtk::ParallelGridCreator<ALUGridType<2>>>(mpi, test, "ALUGridType<2, Parallel>", false);

  reader_writer_test<ALUGridType<3>, Vtk::SerialGridCreator<ALUGridType<3>>>(mpi, test, "ALUGridType<3>");
  #if 0
  // Currently the 2.7 branch is not working, due to a new bisection compatibility check in 3d
  reader_writer_test<ALUGridType<3>, Vtk::ParallelGridCreator<ALUGridType<3>>>(mpi, test, "ALUGridType<3, Parallel>", false);
  #endif
#endif

  return test.exit();
}
