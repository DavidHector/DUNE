#ifndef DUNE_NONLINOPT_TEST_OSBORNE1_HH
#define DUNE_NONLINOPT_TEST_OSBORNE1_HH

#include "leastsquares.hh"

/**
 * @brief Test problem 17 from Moré, Garbow and Hillstrom
 */
class Osborne1Problem
: public TestLeastSquaresProblem
{
  public:

    using Real  = typename TestLeastSquaresProblem::Real;
    using Point = typename TestLeastSquaresProblem::Point;

  private:
    std::vector<Real> y = {0.844,0.908,0.932,0.936,0.925,0.908,
    0.881,0.850,0.818,0.784,0.751,0.718,0.685,0.658,0.628,0.603,
    0.580,0.558,0.538,0.522,0.506,0.490,0.478,0.467,0.457,0.448,
    0.438,0.431,0.424,0.420,0.414,0.411,0.406};

  public:

    Osborne1Problem()
      : TestLeastSquaresProblem(5,33,1,1)
    {}

    std::string name() const override
    {
      return "Osborne 1 function";
    }

    Real component(unsigned int i, const Point& x) const override
    {
      const Real t_i = 10.*i;
      return y[i] - (x[0] + x[1]*std::exp(-t_i*x[3]) + x[2]*std::exp(-t_i*x[4]));
    }

    void gradient_comp(unsigned int i, const Point& x,
        Point& gradient) const override
    {
      const Real t_i = 10.*i;
      gradient[0] += -2. * 1. * component(i,x);
      gradient[1] += -2. * std::exp(-t_i*x[3]) * component(i,x);
      gradient[2] += -2. * std::exp(-t_i*x[4]) * component(i,x);
      gradient[3] +=  2. * t_i*x[1]*std::exp(-t_i*x[3]) * component(i,x);
      gradient[4] +=  2. * t_i*x[2]*std::exp(-t_i*x[4]) * component(i,x);
    }

    Point start_point(Real scale = 1.) const override
    {
      return {0.5*scale,1.5*scale,-scale,0.01*scale,0.02*scale};
    }

    Real goal_value(unsigned int i) const override
    {
      return 5.464894697482767e-05;
    }

    Point goal_point(unsigned int i) const override
    {
      // unknown, best found point
      return {0.37541,1.93585,-1.46469,0.0128675,0.0221227};
    }
};

#endif // DUNE_NONLINOPT_TEST_OSBORNE1_HH
