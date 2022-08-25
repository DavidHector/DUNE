#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <iomanip>

#if HAVE_DUNE_COMMON
#include<dune/common/exceptions.hh>
#endif // HAVE_DUNE_COMMON

#include<dune/nonlinopt/nonlinopt.hh>

#include "leastsquares.hh"

std::tuple<unsigned int, unsigned int, unsigned int, unsigned int,
  unsigned int, bool, double, double, double>
solve(const TestLeastSquaresProblem& problem, double scale, unsigned int max_iter)
{
  typename TestLeastSquaresProblem::Point point(problem.start_point(scale));

  typename Dune::NonlinOpt::UnconstrainedOptimization<double,
        Dune::NonlinOpt::VectorClass<double>> solver(max_iter,10,0);

  try
  {
    solver.apply(problem,point);
  }
  catch(...)
  {
    std::cout << "=== linesearch failed!" << std::endl;
  }

  return problem.report(point,solver.get_value());
}

#if HAVE_DUNE_COMMON
std::tuple<unsigned int, unsigned int, unsigned int, unsigned int,
  unsigned int, bool, double, double, double>
solve(const TestLeastSquaresProblem& problem, double scale,
    const Dune::ParameterTree& config)
{
  typename TestLeastSquaresProblem::Point point(problem.start_point(scale));

  typename Dune::NonlinOpt::UnconstrainedOptimization<double,
        Dune::NonlinOpt::VectorClass<double>> solver(config);

  try
  {
    solver.apply(problem,point);
  }
  catch(...)
  {
    std::cout << "=== linesearch failed!" << std::endl;
  }

  return problem.report(point,solver.get_value());
}
#endif // HAVE_DUNE_COMMON
