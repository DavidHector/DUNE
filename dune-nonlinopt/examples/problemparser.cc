#ifdef HAVE_MUPARSER

#include<dune/nonlinopt/nonlinopt.hh>

int main(int argc, char** argv)
{
  if (argc != 3)
  {
    std::cout << "usage: <string containing expression> <string containing initial point>\n"
      << "the initial point has to be specified as numbers separated by whitespace\n\n"
      << "some possible expressions in 2D:\n"
      << "trivial quadratic: \"x^2+y^2\"\n"
      << "slanted quadratic: \"0.1*(x-y)^2 + 10*(x+y)^2\"\n"
      << "Rosenbrock:        \"(1-x)^2 + 100*(y-x^2)^2\"\n"
      << "Beale:             \"(1.5-x+x*y)^2 + (2.25-x+x*y^2)^2 + (2.625-x+x*y^3)^2\"\n"
      << "Booth:             \"(x+2*y-7)^2 + (2*x+y-5)^2\"\n"
      << "Himmelblau:        \"(x^2+y-11)^2 + (x+y^2-7)^2\"\n";
    return 0;
  }

  // create problem definition
  using Problem = Dune::NonlinOpt::ParsedProblem;
  const std::string expression(argv[1]);
  Problem problem(expression);

  // read in initial position
  using Point = typename Problem::Point;
  Point point = problem.zero();
  const std::string initial(argv[2]);
  std::stringstream init_stream(initial);
  for (std::size_t i = 0; i < problem.dim(); i++)
    init_stream >> point[i];

  // display information about initial position
  std::cout << std::scientific;
  std::cout << std::endl << "initial point: (";
  for (std::size_t i = 0; i < problem.dim()-1; i++)
    std::cout << point[i] << ", ";
  std::cout << point[problem.dim()-1] << ")" << std::endl;

  std::cout << "initial value: " << problem.value(point) << std::endl;
  Point gradient;
  problem.gradient(point,gradient);
  std::cout << "initial gradient: (";
  for (std::size_t i = 0; i < problem.dim()-1; i++)
    std::cout << gradient[i] << ", ";
  std::cout << gradient[problem.dim()-1] << ")" << std::endl;
  std::cout << "two-norm: " << std::sqrt(gradient * gradient)
    << " inf-norm: " << gradient.inf_norm() << std::endl << std::endl;

  // create solver and find local minimum
  using Solver = typename Dune::NonlinOpt::UnconstrainedOptimization<double,Point>;
  Solver solver;

  solver.apply(problem,point);

  // display information about final position
  std::cout << std::endl << "final point: (";
  for (std::size_t i = 0; i < problem.dim()-1; i++)
    std::cout << point[i] << ", ";
  std::cout << point[problem.dim()-1] << ")" << std::endl;

  std::cout << "final value: " << problem.value(point) << std::endl;
  problem.gradient(point,gradient);
  std::cout << "final gradient: (";
  for (std::size_t i = 0; i < problem.dim()-1; i++)
    std::cout << gradient[i] << ", ";
  std::cout << gradient[problem.dim()-1] << ")" << std::endl;
  std::cout << "two-norm: " << std::sqrt(gradient * gradient)
    << " inf-norm: " << gradient.inf_norm() << std::endl;
  std::cout << "function evals: " << problem.evals() - 8 << std::endl;
}

#else // HAVE_MUPARSER
#include <iostream>

int main()
{
  std::cout << "this example requires the muParser library\n";
}
#endif // HAVE_MUPARSER
