#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <iomanip>


#include "${problem_lower}.hh"

int main()
{
  std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int,
    unsigned int, bool, double, double, double>> results;

  Dune::ParameterTree config;
  config["optimization.verbosity"]     = "${verbosity}";
  config["optimization.reset_iter"]    = "${max_iter}";
  config["optimization.reset_iter"]    = "${reset_iter}";
  config["optimization.use_scale"]     = "${use_scale}";
  config["optimization.conjugation"]   = "${conjugation}";
  config["optimization.quasi_newton"]  = "${quasi_newton}";
  config["optimization.extrapolation"] = "${extrapolation}";
  config["optimization.linesearch"]    = "${linesearch}";
  config["storage.limit"]              = "${storage_limit}";
  config["linesearch.interpolation"]   = "${interpolation}";

  ${problem}Problem problem;
  for (double scale = 0.1; scale < 999.; scale *= 10.)
    results.push_back(solve(problem,scale,config));

  std::cout << std::scientific << std::setprecision(6);
  std::cout << "${problem} ${solver}" << std::endl;
  std::cout << "results:" << std::endl;
  std::cout << "ITER  f     g     f+g   f+3g   CONV  DESC_2NORM    RESIDUAL      ERROR" << std::endl;
  std::cout << "-----------------------------------------------------------------------------" << std::endl;
  for (const auto& e : results)
  {
    std::cout << std::setw(5) << std::get<0>(e)
      << " " << std::setw(5) << std::get<1>(e)
      << " " << std::setw(5) << std::get<2>(e)
      << " " << std::setw(5) << std::get<3>(e)
      << " " << std::setw(6) << std::get<4>(e)
      << " " << std::setw(4) << std::get<5>(e)
      << " " << std::setw(13) << std::get<6>(e)
      << " " << std::setw(13) << std::get<7>(e)
      << " " << std::setw(13) << std::get<8>(e) << std::endl;
  }

  if (! std::get<5>(results[1]))
    return 1;

  return 0;
}
