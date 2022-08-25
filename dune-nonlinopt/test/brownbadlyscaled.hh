#ifndef DUNE_NONLINOPT_TEST_BROWNBADLY_HH
#define DUNE_NONLINOPT_TEST_BROWNBADLY_HH

#include "leastsquares.hh"

/**
 * @brief Test problem 4 from Moré, Garbow and Hillstrom
 */
class BrownBadlyScaledProblem
: public TestLeastSquaresProblem
{
  public:

    using Real  = typename TestLeastSquaresProblem::Real;
    using Point = typename TestLeastSquaresProblem::Point;

    BrownBadlyScaledProblem()
      : TestLeastSquaresProblem(2,3)
    {}

    std::string name() const override
    {
      return "Brown badly scaled function";
    }

    Real component(unsigned int i, const Point& x) const override
    {
      if (i == 0)
        return x[0] - 1e6;
      else if (i == 1)
        return x[1] - 2e-6;
      else
        return
          x[0]*x[1] -2.;
    }

    void gradient_comp(unsigned int i, const Point& x,
        Point& gradient) const override
    {
      if (i == 0)
        gradient[0] += 2. * (x[0] - 1e6);
      else if (i == 1)
        gradient[1] += 2. * (x[1] - 2e-6);
      else
      {
        gradient[0] += 2. * x[1] * (x[0]*x[1] - 2.);
        gradient[1] += 2. * x[0] * (x[0]*x[1] - 2.);
      }
    }

    Point start_point(Real scale = 1.) const override
    {
      return {scale,scale};
    }

    Real goal_value(unsigned int i) const override
    {
      return 0.;
    }

    Point goal_point(unsigned int i) const override
    {
      return {1e6,2e-6};
    }
};

#endif // DUNE_NONLINOPT_TEST_BROWNBADLY_HH
