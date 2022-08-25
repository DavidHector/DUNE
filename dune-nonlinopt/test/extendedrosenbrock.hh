#ifndef DUNE_NONLINOPT_TEST_EXTENDEDROSENBROCK_HH
#define DUNE_NONLINOPT_TEST_EXTENDEDROSENBROCK_HH

#include "rosenbrock.hh"

/**
 * @brief Test problem 21 from Moré, Garbow and Hillstrom
 */
class ExtendedRosenbrockProblem
: public RosenbrockProblem
{
  public:

    ExtendedRosenbrockProblem()
      : RosenbrockProblem(1.,10.,100)
    {}
};

#endif // DUNE_NONLINOPT_TEST_EXTENDEDROSENBROCK_HH
