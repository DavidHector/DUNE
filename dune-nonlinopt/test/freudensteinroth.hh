#ifndef DUNE_NONLINOPT_TEST_FREUDENSTEINROTH_HH
#define DUNE_NONLINOPT_TEST_FREUDENSTEINROTH_HH

#include "leastsquares.hh"

/**
 * @brief Test problem 2 from Moré, Garbow and Hillstrom
 */
class FreudensteinRothProblem
: public TestLeastSquaresProblem
{
  public:

    using Real  = typename TestLeastSquaresProblem::Real;
    using Point = typename TestLeastSquaresProblem::Point;

    FreudensteinRothProblem()
      : TestLeastSquaresProblem(2,2,2,2)
    {}

    std::string name() const override
    {
      return "Freudenstein and Roth function";
    }

    Real component(unsigned int i, const Point& x) const override
    {
      if (i == 0)
        return -13. + x[0] + ((5. - x[1])*x[1] - 2.)*x[1];
      else
        return -29. + x[0] + ((x[1] + 1.)*x[1] - 14.)*x[1];
    }

    void gradient_comp(unsigned int i, const Point& x,
        Point& gradient) const override
    {
      if (i == 0)
      {
        gradient[0] +=  2. * (-13. + x[0] + ((5. - x[1])*x[1] - 2.)*x[1]);
        gradient[1] += -2. * ((3.*x[1] - 10.)*x[1] + 2.)
          * (-13. + x[0] + ((5. - x[1])*x[1] - 2.)*x[1]);
      }
      else
      {
        gradient[0] += 2. * (-29. + x[0] + ((x[1] + 1.)*x[1] - 14.)*x[1]);
        gradient[1] += 2. * ((3.*x[1] + 2.)*x[1] - 14.)
          * (-29. + x[0] + ((x[1] + 1.)*x[1] - 14.)*x[1]);
      }
    }

    Point start_point(Real scale = 1.) const override
    {
      return {0.5*scale, -2.*scale};
    }

    Real goal_value(unsigned int i) const override
    {
      if (i == 0)
        return 0.;
      else
        return 48.98425367924000;
    }

    Point goal_point(unsigned int i) const override
    {
      if (i == 0)
        return {5.,4.};
      else
        return {11.4128,-0.896805};
    }
};

#endif // DUNE_NONLINOPT_TEST_FREUDENSTEINROTH_HH
