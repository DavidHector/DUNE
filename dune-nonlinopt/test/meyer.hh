#ifndef DUNE_NONLINOPT_TEST_MEYER_HH
#define DUNE_NONLINOPT_TEST_MEYER_HH

#include "leastsquares.hh"

/**
 * @brief Test problem 10 from Moré, Garbow and Hillstrom
 */
class MeyerProblem
: public TestLeastSquaresProblem
{
  public:

    using Real  = typename TestLeastSquaresProblem::Real;
    using Point = typename TestLeastSquaresProblem::Point;

  private:

    std::vector<Real> y = {34780,28610,23650,19630,16370,
    13720,11540,9744,8261,7030,6005,5147,4427,3820,3307,
    2872};

  public:

    MeyerProblem()
      : TestLeastSquaresProblem(3,16,1,1)
    {}

    std::string name() const override
    {
      return "Meyer function";
    }

    Real component(unsigned int i, const Point& x) const override
    {
      const Real t_i = 45. + 5.*(i+1);
      return x[0] * std::exp(x[1]/(t_i + x[2])) - y[i];
    }

    void gradient_comp(unsigned int i, const Point& x,
        Point& gradient) const override
    {
      const Real t_i = 45. + 5.*(i+1);
      gradient[0] +=  2. * std::exp(x[1]/(t_i + x[2])) * component(i,x);
      gradient[1] +=  2. * x[0]/(t_i + x[2])
        * std::exp(x[1]/(t_i + x[2])) * component(i,x);
      gradient[2] += -2. * x[0]*x[1]/std::pow(t_i + x[2],2)
        * std::exp(x[1]/(t_i + x[2])) * component(i,x);
    }

    Point start_point(Real scale = 1.) const override
    {
      return {0.2*scale,4e3*scale,250.*scale};
    }

    Real goal_value(unsigned int i) const override
    {
      return 8.794585517073886e+01;
    }

    Point goal_point(unsigned int i) const override
    {
      // unknown, best found point
      return {0.00560964,6181.35,345.224};
    }
};

#endif // DUNE_NONLINOPT_TEST_MEYER_HH
