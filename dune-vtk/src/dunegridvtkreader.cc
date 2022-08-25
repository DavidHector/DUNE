#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/vtk/vtkreader.hh>
#include <dune/vtk/gridcreators/discontinuousgridcreator.hh>
#include <dune/vtk/utility/errors.hh>

#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#elif HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

using namespace Dune;

template <class GridType>
void test(std::string prefix, GridType const& grid)
{
  { // write conforming binary file
    Dune::VTKWriter vtkWriter{grid.leafGridView(), Dune::VTK::DataMode::conforming, Dune::VTK::Precision::float32};
    auto fn = vtkWriter.write(prefix + "_conforming_binary", Dune::VTK::OutputType::appendedraw);

    // Read conforming binary file
    std::cout << "read '" << fn << "'..." << std::endl;
    auto gridPtr = VtkReader<GridType>::createGridFromFile(fn);

    VTK_ASSERT(gridPtr->size(0) == grid.size(0));
    VTK_ASSERT(gridPtr->size(GridType::dimension) == grid.size(GridType::dimension));
    std::remove(fn.c_str());
  }

  { // write non-conforming ascii file
    Dune::VTKWriter vtkWriter{grid.leafGridView(), Dune::VTK::DataMode::nonconforming, Dune::VTK::Precision::float32};
    auto fn = vtkWriter.write(prefix + "_nonconforming_ascii", Dune::VTK::OutputType::ascii);

    // Read non-conforming binary file
    std::cout << "read '" << fn << "'..." << std::endl;
    using DGC = Dune::Vtk::DiscontinuousGridCreator<GridType>;
    auto gridPtr = VtkReader<GridType, DGC>::createGridFromFile(fn);

    VTK_ASSERT(gridPtr->size(0) == grid.size(0));
    VTK_ASSERT(gridPtr->size(GridType::dimension) == grid.size(GridType::dimension));
    std::remove(fn.c_str());
  }
}

#if HAVE_DUNE_UGGRID
  template <int dim> using CubeGrid = Dune::UGGrid<dim>;
  template <int dim> using SimplexGrid = Dune::UGGrid<dim>;
#elif HAVE_DUNE_ALUGRID
  template <int dim> using CubeGrid = Dune::ALUGrid<dim,dim,Dune::cube,Dune::nonconforming>;
  template <int dim> using SimplexGrid = Dune::ALUGrid<dim,dim,Dune::simplex,Dune::conforming>;
#endif

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  std::array<unsigned int,2> subdivision2{2u,2u};
  std::array<unsigned int,3> subdivision3{2u,2u,2u};

#if HAVE_DUNE_UGGRID || HAVE_DUNE_ALUGRID
  {
    auto grid = Dune::StructuredGridFactory<CubeGrid<2>>::createCubeGrid({0.0, 0.0}, {1.0, 1.0}, subdivision2);
    test("dunegridvtkreader_cube_2d", *grid);
  }

  {
    auto grid = Dune::StructuredGridFactory<SimplexGrid<2>>::createSimplexGrid({0.0, 0.0}, {1.0, 1.0}, subdivision2);
    test("dunegridvtkreader_simplex_2d", *grid);
  }

  {
    auto grid = Dune::StructuredGridFactory<CubeGrid<3>>::createCubeGrid({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, subdivision3);
    test("dunegridvtkreader_cube_3d", *grid);
  }

  {
    auto grid = Dune::StructuredGridFactory<SimplexGrid<3>>::createSimplexGrid({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, subdivision3);
    test("dunegridvtkreader_simplex_3d", *grid);
  }
#endif
}
