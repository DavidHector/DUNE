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

#include <dune/grid/yaspgrid.hh>

#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>
#include <dune/vtk/legacyvtkfunction.hh>

using namespace Dune;

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  const int dim = 3;
  using GridType = YaspGrid<dim>;
  FieldVector<double,dim> upperRight; upperRight = 1.0;
  auto numElements = filledArray<dim,int>(8);
  GridType grid(upperRight,numElements,0,0);

  using GridView = typename GridType::LeafGridView;
  GridView gridView = grid.leafGridView();

  std::vector<double> p1function(gridView.size(dim), 1.0);
  using P1Function = P1VTKFunction<GridView,std::vector<double>>;
  std::shared_ptr<VTKFunction<GridView> const> p1FctWrapped(new P1Function(gridView, p1function, "p1"));

  using Writer = VtkUnstructuredGridWriter<GridView>;
  Writer vtkWriter(gridView, Vtk::FormatTypes::ASCII);
  vtkWriter.addPointData(p1FctWrapped);
  vtkWriter.write("legacyvtkwriter_ascii_float32.vtu");
}
