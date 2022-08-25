// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <cstring>
#include <iostream>
#include <set>
#include <vector>

#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/filledarray.hh>
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
#include <dune/vtk/gridcreators/continuousgridcreator.hh>

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

using TestCases = std::set<std::tuple<std::string,Vtk::FormatTypes,Vtk::DataTypes>>;
static TestCases test_cases = {
  {"ascii32", Vtk::FormatTypes::ASCII, Vtk::DataTypes::FLOAT32},
  {"bin32", Vtk::FormatTypes::BINARY, Vtk::DataTypes::FLOAT32},
  {"zlib32", Vtk::FormatTypes::COMPRESSED, Vtk::DataTypes::FLOAT32},
  {"ascii64", Vtk::FormatTypes::ASCII, Vtk::DataTypes::FLOAT64},
  {"bin64", Vtk::FormatTypes::BINARY, Vtk::DataTypes::FLOAT64},
  {"zlib64", Vtk::FormatTypes::COMPRESSED, Vtk::DataTypes::FLOAT64},
};

template <class Test>
void compare (Test& test, Vtk::Path const& dir, Vtk::Path const& name)
{
  test.check(compare_files(dir.string() + '/' + name.string() + ".vtu",
                           dir.string() + '/' + name.string() + "_2.vtu"));
}

template <class GridView>
void writer_test (GridView const& gridView)
{
  for (auto const& test_case : test_cases) {
    VtkUnstructuredGridWriter<GridView> vtkWriter(gridView, std::get<1>(test_case), std::get<2>(test_case));
    vtkWriter.write("reader_writer_test_" + std::get<0>(test_case) + ".vtu");
  }
}

template <class Grid, class Test>
void reader_test (MPIHelper& mpi, Test& test)
{
  std::string ext = ".vtu";
  for (auto const& test_case : test_cases) {
    std::vector<std::string> pieces1, pieces2;

    {
      GridFactory<Grid> factory;
      VtkReader<Grid> reader{factory};
      reader.read("reader_writer_test_" + std::get<0>(test_case) + ext);

      std::unique_ptr<Grid> grid{factory.createGrid()};
      pieces1 = reader.pieces();

      VtkUnstructuredGridWriter<typename Grid::LeafGridView> vtkWriter(grid->leafGridView(),
        std::get<1>(test_case), std::get<2>(test_case));
      vtkWriter.write("reader_writer_test_" + std::get<0>(test_case) + "_2.vtu");
    }

    {
      GridFactory<Grid> factory2;
      VtkReader<Grid> reader2{factory2};
      reader2.read("reader_writer_test_" + std::get<0>(test_case) + "_2" + ext, false);
      pieces2 = reader2.pieces();
    }

    test.check(pieces1.size() == pieces2.size(), "pieces1.size == pieces2.size");
    for (std::size_t i = 0; i < pieces1.size(); ++i)
      test.check(compare_files(pieces1[i], pieces2[i]));
  }
}


template <int I>
using int_ = std::integral_constant<int,I>;

int main (int argc, char** argv)
{
  auto& mpi = Dune::MPIHelper::instance(argc, argv);
  if (mpi.size() > 1) {
    std::cout << "Parallel VtkReader not yet supported\n";
    return 0;
  }

  TestSuite test{};

#if HAVE_DUNE_UGGRID
  // Test VtkWriter for UGGrid
  Hybrid::forEach(std::make_tuple(int_<2>{}, int_<3>{}), [&test,&mpi](auto dim)
  {
    using GridType = UGGrid<dim.value>;
    {
      FieldVector<double,dim.value> lowerLeft; lowerLeft = 0.0;
      FieldVector<double,dim.value> upperRight; upperRight = 1.0;
      auto numElements = filledArray<dim.value,unsigned int>(4);
      auto gridPtr = StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft, upperRight, numElements);
      gridPtr->loadBalance();

      writer_test(gridPtr->leafGridView());
    }

    reader_test<GridType>(mpi,test);
  });
#endif

#if 0 && HAVE_DUNE_ALUGRID
  // Test VtkWriter for ALUGrid. Currently the 2.7 branch is not working.
  Hybrid::forEach(std::make_tuple(int_<2>{}, int_<3>{}), [&test,&mpi](auto dim)
  {
    using GridType = Dune::ALUGrid<dim.value, dim.value, Dune::simplex, Dune::conforming>;
    {
      FieldVector<double,dim.value> lowerLeft; lowerLeft = 0.0;
      FieldVector<double,dim.value> upperRight; upperRight = 1.0;
      auto numElements = filledArray<dim.value,unsigned int>(4);
      auto gridPtr = StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft, upperRight, numElements);
      gridPtr->loadBalance();

      writer_test(gridPtr->leafGridView());
    }

    reader_test<GridType>(mpi,test);
  });
#endif

  // Test VtkWriter for YaspGrid
  Hybrid::forEach(std::make_tuple(int_<1>{}, int_<2>{}, int_<3>{}), [](auto dim)
  {
    using GridType = YaspGrid<dim.value>;
    FieldVector<double,dim.value> upperRight; upperRight = 1.0;
    auto numElements = filledArray<dim.value,int>(8);
    GridType grid(upperRight, numElements);
    writer_test(grid.leafGridView());
  });

  return test.exit();
}
