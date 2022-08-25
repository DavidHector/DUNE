#ifndef DUNE_NONLINOPT_TEST_POWELLBADLY_HH
#define DUNE_NONLINOPT_TEST_POWELLBADLY_HH

#include "leastsquares.hh"

/**
 * @brief Test problem 3 from Moré, Garbow and Hillstrom
 */
class PowellBadlyScaledProblem
: public TestLeastSquaresProblem
{
  public:

    using Real  = typename TestLeastSquaresProblem::Real;
    using Point = typename TestLeastSquaresProblem::Point;

    PowellBadlyScaledProblem()
      : TestLeastSquaresProblem(2,2,1,2)
    {}

    std::string name() const override
    {
      return "Powell badly scaled function";
    }

    Real component(unsigned int i, const Point& x) const override
    {
      if (i == 0)
        return 1e4 * x[0]*x[1] - 1.;
      else
        return std::exp(-x[0]) + std::exp(-x[1]) - 1.0001;
    }

    void gradient_comp(unsigned int i, const Point& x,
        Point& gradient) const override
    {
      if (i == 0)
      {
        gradient[0] += 2. * 1e4 * x[1] * (1e4 * x[0]*x[1] - 1.);
        gradient[1] += 2. * 1e4 * x[0] * (1e4 * x[0]*x[1] - 1.);
      }
      else
      {
        gradient[0] += -2. * std::exp(-x[0]) * (std::exp(-x[0]) + std::exp(-x[1]) - 1.0001);
        gradient[1] += -2. * std::exp(-x[1]) * (std::exp(-x[0]) + std::exp(-x[1]) - 1.0001);
      }
    }

    Point start_point(Real scale = 1.) const override
    {
      return {0.,scale};
    }

    Real goal_value(unsigned int i) const override
    {
      return 0.;
    }

    Point goal_point(unsigned int i) const override
    {
      if (i == 0)
        return {1.09816e-05,9.10615};
      else
        return {1e-5,10.};
    }
};

#endif // DUNE_NONLINOPT_TEST_POWELLBADLY_HH
