#ifndef DUNE_NONLINOPT_TEST_GAUSSIAN_HH
#define DUNE_NONLINOPT_TEST_GAUSSIAN_HH

#include "leastsquares.hh"

/**
 * @brief Test problem 9 from Moré, Garbow and Hillstrom
 */
class GaussianProblem
: public TestLeastSquaresProblem
{
  public:

    using Real  = typename TestLeastSquaresProblem::Real;
    using Point = typename TestLeastSquaresProblem::Point;

  private:

    std::vector<Real> y = {0.0009,0.0044,0.0175,0.0540,
      0.1295,0.2420,0.3521,0.3989,0.3521,0.2420,0.1295,
      0.0540,0.0175,0.0044,0.0009};

  public:

    GaussianProblem()
      : TestLeastSquaresProblem(3,15,1,1)
    {}

    std::string name() const override
    {
      return "Gaussian function";
    }

    Real component(unsigned int i, const Point& x) const override
    {
      return x[0] * std::exp(-x[1]*std::pow((8.-(i+1))/2.-x[2],2)/2) - y[i];
    }

    void gradient_comp(unsigned int i, const Point& x,
        Point& gradient) const override
    {
      gradient[0] +=  2. * std::exp(-x[1]*std::pow((8.-(i+1))/2.-x[2],2)/2) * component(i,x);
      gradient[1] += -2. * x[0]*std::pow((8.-(i+1))/2.-x[2],2)/2
        * std::exp(-x[1]*std::pow((8.-(i+1))/2.-x[2],2)/2) * component(i,x);
      gradient[2] +=  2. * x[0]*x[1]*((8.-(i+1))/2.-x[2])
        * std::exp(-x[1]*std::pow((8.-(i+1))/2.-x[2],2)/2) * component(i,x);
    }

    Point start_point(Real scale = 1.) const override
    {
      return {0.4*scale,scale,0.};
    }

    Real goal_value(unsigned int i) const override
    {
      return 1.127932769699661e-08;
    }

    Point goal_point(unsigned int i) const override
    {
      // unknown, best found point
      return {0.398956,1.00002,0.};
    }
};

#endif // DUNE_NONLINOPT_TEST_GAUSSIAN_HH
