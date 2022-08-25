#ifndef DUNE_NONLINOPT_TEST_OSBORNE2_HH
#define DUNE_NONLINOPT_TEST_OSBORNE2_HH

#include "leastsquares.hh"

/**
 * @brief Test problem 19 from Moré, Garbow and Hillstrom
 */
class Osborne2Problem
: public TestLeastSquaresProblem
{
  public:

    using Real  = typename TestLeastSquaresProblem::Real;
    using Point = typename TestLeastSquaresProblem::Point;

  private:
    std::vector<Real> y = {1.366,1.191,1.112,1.013,0.991,0.885,
    0.831,0.847,0.786,0.725,0.746,0.679,0.608,0.655,0.616,0.606,
    0.602,0.626,0.651,0.724,0.649,0.649,0.694,0.644,0.624,0.661,
    0.612,0.558,0.533,0.495,0.500,0.423,0.395,0.375,0.372,0.391,
    0.396,0.405,0.428,0.429,0.523,0.562,0.607,0.653,0.672,0.708,
    0.633,0.668,0.645,0.632,0.591,0.559,0.597,0.625,0.739,0.710,
    0.729,0.720,0.636,0.581,0.428,0.292,0.162,0.098,0.054};

  public:

    Osborne2Problem()
      : TestLeastSquaresProblem(11,65,1,1)
    {}

    std::string name() const override
    {
      return "Osborne 2 function";
    }

    Real component(unsigned int i, const Point& x) const override
    {
      const Real t_i = i/10.;
      return y[i] - (x[0]*std::exp(-t_i*x[4]) + x[1]*std::exp(-std::pow(t_i-x[8],2)*x[5])
          + x[2]*std::exp(-std::pow(t_i-x[9],2)*x[6]) + x[3]*std::exp(-std::pow(t_i-x[10],2)*x[7]));
    }

    void gradient_comp(unsigned int i, const Point& x,
        Point& gradient) const override
    {
      const Real t_i = i/10.;
      gradient[0]  += -2. * std::exp(-t_i*x[4]) * component(i,x);
      gradient[1]  += -2. * std::exp(-std::pow(t_i-x[8],2)*x[5]) * component(i,x);
      gradient[2]  += -2. * std::exp(-std::pow(t_i-x[9],2)*x[6]) * component(i,x);
      gradient[3]  += -2. * std::exp(-std::pow(t_i-x[10],2)*x[7]) * component(i,x);
      gradient[4]  +=  2. * t_i*x[0]*std::exp(-t_i*x[4]) * component(i,x);
      gradient[5]  +=  2. * std::pow(t_i-x[8],2)*x[1]*std::exp(-std::pow(t_i-x[8],2)*x[5]) * component(i,x);
      gradient[6]  +=  2. * std::pow(t_i-x[9],2)*x[2]*std::exp(-std::pow(t_i-x[9],2)*x[6]) * component(i,x);
      gradient[7]  +=  2. * std::pow(t_i-x[10],2)*x[3]*std::exp(-std::pow(t_i-x[10],2)*x[7]) * component(i,x);
      gradient[8]  += -2. * x[1]*x[5]*2.*(t_i-x[8])*std::exp(-std::pow(t_i-x[8],2)*x[5]) * component(i,x);
      gradient[9]  += -2. * x[2]*x[6]*2.*(t_i-x[9])*std::exp(-std::pow(t_i-x[9],2)*x[6]) * component(i,x);
      gradient[10] += -2. * x[3]*x[7]*2.*(t_i-x[10])*std::exp(-std::pow(t_i-x[10],2)*x[7]) * component(i,x);
    }

    Point start_point(Real scale = 1.) const override
    {
      return {1.3*scale,0.65*scale,0.65*scale,0.7*scale,0.6*scale,
      3.*scale,5.*scale,7.*scale,2.*scale,4.5*scale,5.5*scale};
    }

    Real goal_value(unsigned int i) const override
    {
      return 4.013773629354807e-2;
    }

    Point goal_point(unsigned int i) const override
    {
      // unknown, best found point
      return {1.30998,0.431554,0.633662,0.599431,0.754183,0.904289,1.36581,4.8237,2.39868,4.56887,5.67534};
    }
};

#endif // DUNE_NONLINOPT_TEST_OSBORNE2_HH
