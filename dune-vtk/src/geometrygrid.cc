// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <vector>

#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/geometrygrid.hh>

#include <dune/vtk/vtkwriter.hh>

using namespace Dune;
using namespace Dune::Functions;

class TorusMapper
    : public AnalyticalCoordFunction<double, 2, 3, TorusMapper>
{
  using Super = AnalyticalCoordFunction<double, 2, 3, TorusMapper>;

public:
  using DomainVector = Super::DomainVector;
  using RangeVector = Super::RangeVector;

  TorusMapper(double R, double r)
    : R_(R), r_(r) {}

  void evaluate(DomainVector const& x, RangeVector& y) const
  {
    y[0] = (R_ + r_*std::cos(x[0])) * std::cos(x[1]);
    y[1] = (R_ + r_*std::cos(x[0])) * std::sin(x[1]);
    y[2] = r_*std::sin(x[0]);
  }

private:
  double R_, r_;
};

template <class GridView>
void write (std::string prefix, GridView const& gridView)
{
  FieldVector<double,GridView::dimensionworld> c;
  if (GridView::dimensionworld > 0) c[0] = 11.0;
  if (GridView::dimensionworld > 1) c[1] = 7.0;
  if (GridView::dimensionworld > 2) c[2] = 3.0;

  auto p1Analytic = makeAnalyticGridViewFunction([&c](auto const& x) { return c.dot(x); }, gridView);

  VtkWriter<GridView> vtkWriter(gridView, Vtk::FormatTypes::ASCII, Vtk::DataTypes::FLOAT32);
  vtkWriter.addPointData(p1Analytic, "q1");
  vtkWriter.addCellData(p1Analytic, "q0");
  vtkWriter.write(prefix + "_ascii.vtu");
}

int main (int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  using HostGrid = YaspGrid<2>;
  FieldVector<double,2> bbox = {2.0*M_PI, 2.0*M_PI};
  std::array<int,2> num = {4, 12};
  HostGrid hostGrid{bbox, num, 0, 0}; //, std::bitset<2>{"11"}};

  // grid build up of mapped coordinates
  double R = 1.0, r = 0.25;
  TorusMapper mapper{R,r};
  using Grid = GeometryGrid<HostGrid,TorusMapper>;
  Grid grid{hostGrid, mapper};

  write("geometrygrid_torus", grid.leafGridView());
}