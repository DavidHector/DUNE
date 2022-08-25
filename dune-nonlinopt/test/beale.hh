#ifndef DUNE_NONLINOPT_TEST_BEALE_HH
#define DUNE_NONLINOPT_TEST_BEALE_HH

#include "leastsquares.hh"

/**
 * @brief Test problem 5 from Moré, Garbow and Hillstrom
 */
class BealeProblem
: public TestLeastSquaresProblem
{
  public:

    using Real  = typename TestLeastSquaresProblem::Real;
    using Point = typename TestLeastSquaresProblem::Point;

    BealeProblem()
      : TestLeastSquaresProblem(2,3)
    {}

    std::string name() const override
    {
      return "Beale function";
    }

    Real component(unsigned int i, const Point& x) const override
    {
      if (i == 0)
        return 1.5   - x[0] * (1. - x[1]);
      else if (i == 1)
        return 2.25  - x[0] * (1. - std::pow(x[1],2));
      else
        return 2.625 - x[0] * (1. - std::pow(x[1],3));
    }

    void gradient_comp(unsigned int i, const Point& x,
        Point& gradient) const override
    {
      if (i == 0)
      {
        gradient[0] += -2. * (1. - x[1])
          * (1.5 - x[0] * (1. - x[1]));
        gradient[1] +=  2. * x[0]
          * (1.5 - x[0] * (1. - x[1]));
      }
      else if (i == 1)
      {
        gradient[0] += -2. * (1. - std::pow(x[1],2))
          * (2.25  - x[0] * (1. - std::pow(x[1],2)));
        gradient[1] +=  2. * 2.*x[1]*x[0]
          * (2.25  - x[0] * (1. - std::pow(x[1],2)));
      }
      else
      {
        gradient[0] += -2. * (1. - std::pow(x[1],3))
          * (2.625 - x[0] * (1. - std::pow(x[1],3)));
        gradient[1] +=  2. * 3.*std::pow(x[1],2)*x[0]
          * (2.625 - x[0] * (1. - std::pow(x[1],3)));
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
      return {3.,0.5};
    }
};

#endif // DUNE_NONLINOPT_TEST_BEALE_HH
