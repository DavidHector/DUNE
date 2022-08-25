#ifndef DUNE_NONLINOPT_TEST_WOOD_HH
#define DUNE_NONLINOPT_TEST_WOOD_HH

#include "leastsquares.hh"

/**
 * @brief Test problem 14 from Moré, Garbow and Hillstrom
 */
class WoodProblem
: public TestLeastSquaresProblem
{
  public:

    using Real  = typename TestLeastSquaresProblem::Real;
    using Point = typename TestLeastSquaresProblem::Point;

    WoodProblem()
      : TestLeastSquaresProblem(4,6)
    {}

    std::string name() const override
    {
      return "Wood function";
    }

    Real component(unsigned int i, const Point& x) const override
    {
      if (i == 0)
        return 10.*(x[1] - std::pow(x[0],2));
      else if (i == 1)
        return 1 - x[0];
      else if (i == 2)
        return std::sqrt(90.) * (x[3] - std::pow(x[2],2));
      else if (i == 3)
        return 1 - x[2];
      else if (i == 4)
        return std::sqrt(10.) * (x[1] + x[3] - 2.);
      else
        return 1./std::sqrt(10.) * (x[1] - x[3]);
    }

    void gradient_comp(unsigned int i, const Point& x,
        Point& gradient) const override
    {
      if (i == 0)
      {
        gradient[0] += -2. * 20.*x[0] * component(i,x);
        gradient[1] +=  2. * 10. * component(i,x);
      }
      else if (i == 1)
        gradient[0] += -2. * component(i,x);
      else if (i == 2)
      {
        gradient[2] += -2. * std::sqrt(90.)*2.*x[2] * component(i,x);
        gradient[3] +=  2. * std::sqrt(90.) * component(i,x);
      }
      else if (i == 3)
        gradient[2] += -2. * component(i,x);
      else if (i == 4)
      {
        gradient[1] += 2. * std::sqrt(10.) * component(i,x);
        gradient[3] += 2. * std::sqrt(10.) * component(i,x);
      }
      else
      {
        gradient[1] +=  2. * 1./std::sqrt(10.) * component(i,x);
        gradient[3] += -2. * 1./std::sqrt(10.) * component(i,x);
      }
    }

    Point start_point(Real scale = 1.) const override
    {
      return {-3.*scale,-scale,-3.*scale,-scale};
    }

    Real goal_value(unsigned int i) const override
    {
      return 0.;
    }

    Point goal_point(unsigned int i) const override
    {
      return {1.,1.,1.,1.};
    }
};

#endif // DUNE_NONLINOPT_TEST_WOOD_HH
