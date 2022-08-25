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
#include <dune/common/test/testsuite.hh>

#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/vtk/vtkreader.hh>
#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>

using namespace Dune;

// see https://stackoverflow.com/questions/6163611/compare-two-files
bool compare_files (std::string const& fn1, std::string const& fn2)
{
  std::ifstream in1(fn1, std::ios::binary);
  std::ifstream in2(fn2, std::ios::binary);

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

template <class GridView>
void writer_test (GridView const& gridView)
{
  for (auto const& test_case : test_cases) {
    VtkUnstructuredGridWriter<GridView> vtkWriter(gridView, std::get<1>(test_case), std::get<2>(test_case));
    vtkWriter.write("mixed_element_test_" + std::get<0>(test_case) + ".vtu");
  }
}

template <class Grid, class Test>
void reader_test (Test& test)
{
  for (auto const& test_case : test_cases) {
    auto grid = VtkReader<Grid>::createGridFromFile("mixed_element_test_" + std::get<0>(test_case) + ".vtu");
    VtkUnstructuredGridWriter<typename Grid::LeafGridView> vtkWriter(grid->leafGridView(),
      std::get<1>(test_case), std::get<2>(test_case));
    vtkWriter.write("mixed_element_test_" + std::get<0>(test_case) + "_2.vtu");
    test.check(compare_files("mixed_element_test_" + std::get<0>(test_case) + ".vtu",
                             "mixed_element_test_" + std::get<0>(test_case) + "_2.vtu"), std::get<0>(test_case));
  }
}


template <int I>
using int_ = std::integral_constant<int,I>;

int main (int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  TestSuite test{};

#ifdef HAVE_UG
  // Test VtkWriter for UGGrid
  using GridType = UGGrid<2>;
  GridFactory<GridType> factory;

  using X = FieldVector<double,2>;
  using E = std::vector<unsigned int>;
  factory.insertVertex(X{0.0, 0.0}); // 0
  factory.insertVertex(X{1.0, 0.0}); // 1
  factory.insertVertex(X{1.0, 1.0}); // 2
  factory.insertVertex(X{0.0, 1.0}); // 3
  factory.insertVertex(X{1.5, 0.5}); // 4

  factory.insertElement(GeometryTypes::quadrilateral, E{0,1,3,2}); // quadrilateral
  factory.insertElement(GeometryTypes::triangle, E{1,4,2}); // triangle

  {
    std::unique_ptr<GridType> gridPtr{ factory.createGrid() };
    writer_test(gridPtr->leafGridView());
  }
  reader_test<GridType>(test);
#endif

  return test.exit();
}
