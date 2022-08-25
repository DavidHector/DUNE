#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <iomanip>

#if HAVE_DUNE_COMMON
#include<dune/common/parametertree.hh>
#include<dune/common/parametertreeparser.hh>
#endif // HAVE_DUNE_COMMON

#include<test/leastsquares.hh>
#include<test/rosenbrock.hh>

/**
 * @brief Print optimization summary
 */
void print_results(const RosenbrockProblem& problem,
    typename RosenbrockProblem::Point& point,
    typename RosenbrockProblem::Real val)
{
  const auto results = problem.report(point,val);
  std::cout << "iterations: " << std::get<0>(results)
    << " f_eval: " << std::get<1>(results)
    << " g_eval: " << std::get<2>(results)
    << " f+g: " << std::get<3>(results)
    << " f+3g: " << std::get<4>(results)
    << std::endl;
  std::cout << "converged: " << std::get<5>(results)
    << " grad_norm: " << std::get<6>(results)
    << " residual: " << std::get<7>(results)
    << " error: " << std::get<8>(results)
    << std::endl;
}

/**
 * @brief Example application 1
 *
 * Creates an L-BFGS solver and solves the problem.
 */
void driver1(const RosenbrockProblem& problem,
    typename RosenbrockProblem::Point& point)
{
  std::cout << "using the default solver configuration" << std::endl;

  // create default solver (L-BFGS)
  using Solver = Dune::NonlinOpt::UnconstrainedOptimization<double,
        Dune::NonlinOpt::VectorClass<double>>;
  Solver solver;

  // print configuration
  solver.report();

  // solve optimization problem
  solver.apply(problem,point);

  // print results
  print_results(problem,point,solver.get_value());
}

/**
 * @brief Example application 2
 *
 * Like driver1, but switches over to nonlinear CG
 */
void driver2(const RosenbrockProblem& problem,
    typename RosenbrockProblem::Point& point)
{
  std::cout << "configuring a conjugate gradients method" << std::endl;

  using Solver = Dune::NonlinOpt::UnconstrainedOptimization<double,
        Dune::NonlinOpt::VectorClass<double>>;
  Solver solver(1000,1,2);

  // configure solver as nonlinear CG
  solver.set_cg();

  // print configuration
  solver.report();

  // solve optimization problem
  solver.apply(problem,point);

  // print results
  print_results(problem,point,solver.get_value());
}

/**
 * @brief Example application 3
 *
 * Demonstrates how different parts of the optimization
 * routine can be exchanged, like the specific choice of
 * conjugate gradients, the line search, or the termination
 * criterion.
 */
void driver3(const RosenbrockProblem& problem,
    typename RosenbrockProblem::Point& point)
{
  std::cout << "configuring a custom method" << std::endl;

  using Solver = Dune::NonlinOpt::UnconstrainedOptimization<double,
        Dune::NonlinOpt::VectorClass<double>>;
  Solver solver(1000,1,2);
  
  // configure solver as nonlinear CG
  solver.set_cg();

  // use Fletcher-Reeves instead of Hager-Zhang
  solver.set_conjugation<Dune::NonlinOpt::FletcherReevesBeta>();

  // choose a different line search
  solver.set_linesearch<Dune::NonlinOpt::CubicLinesearch>();

  // choose appropriate line search criterion
  solver.set_linesearchcriterion<Dune::NonlinOpt::WolfeConditions>();

  // use Euclidean gradient norm as criterion,
  // disable relative tolerance by setting it zero
  solver.set_stoppingcriterion<Dune::NonlinOpt::GradientTwoNormCriterion>(1e-6,0.);

  // print configuration
  solver.report();

  // solve optimization problem
  solver.apply(problem,point);

  // print results
  print_results(problem,point,solver.get_value());
}

/**
 * @brief Example application 4
 *
 * Demonstrates how the method can be configured using
 * Dune::ParameterTree, with the config created inside
 * the program. This function requires dune-common.
 */
void driver4(const RosenbrockProblem& problem,
    typename RosenbrockProblem::Point& point)
{
#if HAVE_DUNE_COMMON
  std::cout << "creating Dune::ParameterTree configuration inside the program" << std::endl;

  // define configuration
  // see dune-common documentation for additional info
  Dune::ParameterTree config;
  config["optimization.verbosity"]  = "2";
  config["optimization.reset_iter"] = "20";
  config["optimization.linesearch"] = "cubic";
  config["storage.limit"]           = "5";
  config["linesearch.verbosity"]    = "2";
  config["linesearch.criterion"]    = "wolfe";

  // create solver with configuration
  using Solver = Dune::NonlinOpt::UnconstrainedOptimization<double,
        Dune::NonlinOpt::VectorClass<double>>;
  Solver solver(config);

  // print configuration
  solver.report();

  // solve optimization problem
  solver.apply(problem,point);

  // print results
  print_results(problem,point,solver.get_value());
#else
  std::cout << "this choice requires dune-common." << std::endl;
#endif // HAVE_DUNE_COMMON
}

/**
 * @brief Example application 5
 *
 * Demonstrates how the method can be configured using
 * Dune::ParameterTree, with the config read from file.
 * This function requires dune-common.
 */
void driver5(const RosenbrockProblem& problem,
    typename RosenbrockProblem::Point& point)
{
#if HAVE_DUNE_COMMON
  std::cout << "reading Dune::ParameterTree configuration from file" << std::endl;

  // define configuration
  // see dune-common documentation for additional info
  Dune::ParameterTree config;
  Dune::ParameterTreeParser parser;
  parser.readINITree("rosenbrock.ini",config);

  // create solver with configuration
  using Solver = Dune::NonlinOpt::UnconstrainedOptimization<double,
        Dune::NonlinOpt::VectorClass<double>>;
  Solver solver(config);

  // print configuration
  solver.report();

  // solve optimization problem
  solver.apply(problem,point);

  // print results
  print_results(problem,point,solver.get_value());
#else
  std::cout << "this choice requires dune-common." << std::endl;
#endif // HAVE_DUNE_COMMON
}

int main(int argc, char** argv)
{
  if (argc != 2)
  {
    std::cout << "usage: <number between 1 and 5>\n\n"
      << "1) apply default method: L-BFGS\n"
      << "2) use nonlinear conjugate gradients\n"
      << "3) use CG with custom linesearch and criterion\n"
      << "4) use configuration from program (requires dune-common)\n"
      << "5) use configuration from file (requires dune-common)\n";
    return 0;
  }

  RosenbrockProblem problem;
  auto point = problem.start_point(1.);

  int i = atoi(argv[1]);
  if (i == 1)
    driver1(problem,point);
  else if (i == 2)
    driver2(problem,point);
  else if (i == 3)
    driver3(problem,point);
  else if (i == 4)
    driver4(problem,point);
  else if (i == 5)
    driver5(problem,point);
  else
    std::cout << "choice not understood, must be between 1 and 5";
}
