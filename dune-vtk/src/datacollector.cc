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
#include <dune/common/version.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>

#include <dune/vtk/datacollectors/continuousdatacollector.hh>
#include <dune/vtk/datacollectors/discontinuousdatacollector.hh>
#include <dune/vtk/datacollectors/quadraticdatacollector.hh>
#include <dune/vtk/datacollectors/lagrangedatacollector.hh>

using namespace Dune;
using namespace Dune::Functions;

template <class DataCollector, class GridView, class Fct1, class Fct2>
void write_dc (std::string prefix, GridView const& gridView, Fct1 const& fct1, Fct2 const& fct2)
{
  VtkUnstructuredGridWriter<GridView, DataCollector> vtkWriter(gridView, Vtk::FormatTypes::ASCII, Vtk::DataTypes::FLOAT32);
  vtkWriter.addPointData(fct1, "p1");
  vtkWriter.addCellData(fct1, "p0");
  vtkWriter.addPointData(fct2, "q1");
  vtkWriter.addCellData(fct2, "q0");

  vtkWriter.write(prefix + "_" + std::to_string(GridView::dimensionworld) + "d_ascii.vtu");
}

template <class GridView>
void write (std::string prefix, GridView const& gridView)
{
  std::cout << prefix << "..." << std::endl;
  using namespace BasisFactory;
  auto basis = makeBasis(gridView, lagrange<1>());

  FieldVector<double,GridView::dimensionworld> c;
  if (GridView::dimensionworld > 0) c[0] = 11.0;
  if (GridView::dimensionworld > 1) c[1] = 7.0;
  if (GridView::dimensionworld > 2) c[2] = 3.0;

  std::vector<double> vec(basis.dimension());
  interpolate(basis, vec, [&c](auto const& x) { return c.dot(x); });

  // write discrete global-basis function
  auto p1Interpol = makeDiscreteGlobalBasisFunction<double>(basis, vec);

  // write analytic function
  auto p1Analytic = makeAnalyticGridViewFunction([&c](auto const& x) { return c.dot(x); }, gridView);

  write_dc<Vtk::ContinuousDataCollector<GridView>>(prefix + "_continuous", gridView, p1Interpol, p1Analytic);
  write_dc<Vtk::DiscontinuousDataCollector<GridView>>(prefix + "_discontinuous", gridView, p1Interpol, p1Analytic);
  write_dc<Vtk::QuadraticDataCollector<GridView>>(prefix + "_quadratic", gridView, p1Interpol, p1Analytic);

  Hybrid::forEach(StaticIntegralRange<int,7,1>{}, [&](auto p) {
    write_dc<Vtk::LagrangeDataCollector<GridView,p>>(prefix + "_lagrange_p" + std::to_string(p), gridView, p1Interpol, p1Analytic);
  });
}

template <int I>
using int_ = std::integral_constant<int,I>;

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  Hybrid::forEach(StaticIntegralRange<int,4,1>{}, [](auto dim)
  {
    using GridType = YaspGrid<dim.value>;
    FieldVector<double,dim.value> upperRight; upperRight = 1.0;
    auto numElements = filledArray<dim.value,int>(2);
    GridType grid(upperRight, numElements, 0, 0);
    write("datacollector_yasp", grid.leafGridView());
  });

#if HAVE_DUNE_UGGRID
  Hybrid::forEach(StaticIntegralRange<int,4,2>{}, [](auto dim)
  {
    using GridType = UGGrid<dim.value>;
    FieldVector<double,dim.value> lowerLeft; lowerLeft = 0.0;
    FieldVector<double,dim.value> upperRight; upperRight = 1.0;
    auto numElements = filledArray<dim.value,unsigned int>(4);
    auto grid = StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft, upperRight, numElements);
    write("datacollector_ug", grid->leafGridView());
  });
#endif
}
