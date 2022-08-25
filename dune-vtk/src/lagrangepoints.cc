// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <vector>

#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/indices.hh>

#include <dune/vtk/utility/lagrangepoints.hh>
#include <dune/localfunctions/lagrange/interpolation.hh>
#include <dune/localfunctions/lagrange/lagrangebasis.hh>
#include <dune/localfunctions/lagrange/lagrangecoefficients.hh>
#include <dune/localfunctions/utility/localfiniteelement.hh>

using namespace Dune;

template <std::size_t dim>
void write (std::string prefix, index_constant<dim>)
{
  for (int order = 1; order < 6; ++order) {
    std::cout << "order: " << order << std::endl;
    {
      Vtk::LagrangePointSet<double, dim> pointSet(order);
      pointSet.build(GeometryTypes::cube(dim));

      std::size_t i = 0;
      std::cout << "Cube:" << GeometryTypes::cube(dim) << std::endl;
      for (auto const& p : pointSet) {
        std::cout << i++ << ") p = " << p.point() << ", key = " << p.localKey() << std::endl;
      }

      using BasisF = LagrangeBasisFactory<Vtk::LagrangePointSet, dim, double, double>;
      using CoefficientF = LagrangeCoefficientsFactory<Vtk::LagrangePointSet, dim, double>;
      using InterpolationF = LagrangeInterpolationFactory<Vtk::LagrangePointSet, dim, double>;
      GenericLocalFiniteElement<BasisF, CoefficientF, InterpolationF> localFE(GeometryTypes::cube(dim), order);

      [[maybe_unused]] auto const& localBasis = localFE.localBasis();
      [[maybe_unused]] auto const& localCoefficints = localFE.localCoefficients();
      [[maybe_unused]] auto const& localInterpolation = localFE.localInterpolation();
    }

    {
      Vtk::LagrangePointSet<double, dim> pointSet(order);
      pointSet.build(GeometryTypes::simplex(dim));

      std::size_t i = 0;
      std::cout << "Simplex:" << GeometryTypes::simplex(dim) << std::endl;
      for (auto const& p : pointSet) {
        std::cout << i++ << ") p = " << p.point() << ", key = " << p.localKey() << std::endl;
      }
    }
  }
}

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  Hybrid::forEach(std::make_tuple(index_constant<2>{}), [](auto dim)
  {
    write("lagrangepoints", dim);
  });
}
