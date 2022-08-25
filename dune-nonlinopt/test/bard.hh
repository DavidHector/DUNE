#ifndef DUNE_NONLINOPT_TEST_BARD_HH
#define DUNE_NONLINOPT_TEST_BARD_HH

#include "leastsquares.hh"

/**
 * @brief Test problem 8 from Moré, Garbow and Hillstrom
 */
class BardProblem
: public TestLeastSquaresProblem
{
  public:

    using Real  = typename TestLeastSquaresProblem::Real;
    using Point = typename TestLeastSquaresProblem::Point;

  private:

    std::vector<Real> y = {0.14,0.18,0.22,0.25,0.29,0.32,
      0.35,0.39,0.37,0.58,0.73,0.96,1.34,2.10,4.39};

  public:

    BardProblem()
      : TestLeastSquaresProblem(3,15,2,1)
    {}

    std::string name() const override
    {
      return "Bard function";
    }

    Real component(unsigned int i, const Point& x) const override
    {
      Real u = i+1;
      Real v = 16. - (i+1);
      Real w = std::min(u,v);

      return y[i] - (x[0] + u/(v*x[1] + w*x[2]));
    }

    void gradient_comp(unsigned int i, const Point& x,
        Point& gradient) const override
    {
      Real u = i+1;
      Real v = 16. - (i+1);
      Real w = std::min(u,v);

      gradient[0] += -2. * (y[i] - (x[0] + u/(v*x[1] + w*x[2])));
      gradient[1] +=  2. * u*v*std::pow(v*x[1] + w*x[2],-2)
        * (y[i] - (x[0] + u/(v*x[1] + w*x[2])));
      gradient[2] +=  2. * u*w*std::pow(v*x[1] + w*x[2],-2)
        * (y[i] - (x[0] + u/(v*x[1] + w*x[2])));
    }

    Point start_point(Real scale = 1.) const override
    {
      return {scale,scale,scale};
    }

    Real goal_value(unsigned int i) const override
    {
      if (i == 0)
        return 8.214877306580347e-03;
      else
        return 17.4286;
    }

    Point goal_point(unsigned int i) const override
    {
      // unknown, best found point
      return {0.0824105,1.13304,2.3437};
    }
};

#endif // DUNE_NONLINOPT_TEST_BARD_HH
