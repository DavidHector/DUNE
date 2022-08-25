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

#include <dune/polygongrid/grid.hh>
#include <dune/polygongrid/gridfactory.hh>

#include <dune/vtk/legacyvtkfunction.hh>
#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>

using namespace Dune;

using GridType = Dune::PolygonGrid<double>;

std::unique_ptr<GridType> createArbitraryGrid ()
{
  const std::vector< Dune::FieldVector< double, 2 > > vertices
    = { { 0.0, 0.0 }, { 0.5, 0.0 }, { 1.0, 0.0 },
        { 0.0, 0.4 }, { 0.5, 0.2 }, { 0.7, 0.4 }, { 1.0, 0.4 },
        { 0.0, 0.7 }, { 0.7, 0.6 }, { 1.0, 0.6 },
        { 0.0, 1.0 }, { 0.3, 1.0 }, { 1.0, 1.0 } };
  const std::vector< std::vector< unsigned int > > polys
    = { { 0, 1, 4, 3 }, { 1, 2, 6, 5, 4 }, { 3, 4, 5, 8, 11, 7 }, { 5, 6, 9, 8 }, { 7, 11, 10 }, { 8, 9, 12, 11 } };

  Dune::GridFactory<GridType> factory;
  for( const auto &vertex : vertices )
    factory.insertVertex( vertex );
  for( const auto &poly : polys )
    factory.insertElement( Dune::GeometryTypes::none( 2 ), poly );
  return std::unique_ptr<GridType>( factory.createGrid() );
}


int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  auto gridPtr = createArbitraryGrid();

  using GridView = typename GridType::LeafGridView;
  GridView gridView = gridPtr->leafGridView();

  using Writer = VtkUnstructuredGridWriter<GridView>;
  Writer vtkWriter(gridView, Vtk::FormatTypes::ASCII);

  std::vector<double> p1function(gridView.size(GridView::dimension), 1.0);
  using P1Function = P1VTKFunction<GridView,std::vector<double>>;
  std::shared_ptr<VTKFunction<GridView> const> p1FctWrapped(new P1Function(gridView, p1function, "p1"));

  std::vector<double> p0function(gridView.size(0), 1.0);
  using P0Function = P0VTKFunction<GridView,std::vector<double>>;
  std::shared_ptr<VTKFunction<GridView> const> p0FctWrapped(new P0Function(gridView, p0function, "p0"));

  vtkWriter.addPointData(p1FctWrapped);
  vtkWriter.addCellData(p0FctWrapped);

  vtkWriter.write("poly_ascii_float32.vtu");
}
