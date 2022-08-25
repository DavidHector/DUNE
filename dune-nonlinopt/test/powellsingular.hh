#ifndef DUNE_NONLINOPT_TEST_POWELLSINGULAR_HH
#define DUNE_NONLINOPT_TEST_POWELLSINGULAR_HH

#include "leastsquares.hh"

/**
 * @brief Test problem 13 from Moré, Garbow and Hillstrom
 */
class PowellSingularProblem
: public TestLeastSquaresProblem
{
  public:

    using Real  = typename TestLeastSquaresProblem::Real;
    using Point = typename TestLeastSquaresProblem::Point;

    PowellSingularProblem()
      : TestLeastSquaresProblem(4,4)
    {}

    std::string name() const override
    {
      return "Powell singular function";
    }

    Real component(unsigned int i, const Point& x) const override
    {
      if (i == 0)
        return x[0] + 10.*x[1];
      else if (i == 1)
        return std::sqrt(5.)*(x[2] - x[3]);
      else if (i == 2)
        return std::pow(x[1] - 2.*x[2],2);
      else
        return std::sqrt(10.)*std::pow(x[0]-x[3],2);
    }

    void gradient_comp(unsigned int i, const Point& x,
        Point& gradient) const override
    {
      if (i == 0)
      {
        gradient[0] += 2. * component(i,x);
        gradient[1] += 2. * 10. * component(i,x);
      }
      else if (i == 1)
      {
        gradient[2] +=  2. * std::sqrt(5.) * component(i,x);
        gradient[3] += -2. * std::sqrt(5.) * component(i,x);
      }
      else if (i == 2)
      {
        gradient[1] +=  2. * 2.*(x[1] - 2.*x[2]) * component(i,x);
        gradient[2] += -2. * 4.*(x[1] - 2.*x[2]) * component(i,x);
      }
      else
      {
        gradient[0] +=  2. * std::sqrt(10.)*2.*(x[0]-x[3]) * component(i,x);
        gradient[3] += -2. * std::sqrt(10.)*2.*(x[0]-x[3]) * component(i,x);
      }
    }

    Point start_point(Real scale = 1.) const override
    {
      return {3.*scale,-scale,0.,scale};
    }

    Real goal_value(unsigned int i) const override
    {
      return 0.;
    }

    Point goal_point(unsigned int i) const override
    {
      return {0.,0.,0.,0.};
    }
};

#endif // DUNE_NONLINOPT_TEST_POWELLSINGULAR_HH
