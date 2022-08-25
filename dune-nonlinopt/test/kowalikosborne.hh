#ifndef DUNE_NONLINOPT_TEST_KOWALIKOSBORNE_HH
#define DUNE_NONLINOPT_TEST_KOWALIKOSBORNE_HH

#include "leastsquares.hh"

/**
 * @brief Test problem 15 from Moré, Garbow and Hillstrom
 */
class KowalikOsborneProblem
: public TestLeastSquaresProblem
{
  public:

    using Real  = typename TestLeastSquaresProblem::Real;
    using Point = typename TestLeastSquaresProblem::Point;

  private:
    std::vector<Real> y = {0.1957,0.1947,0.1735,0.1600,0.0844,
      0.0627,0.0456,0.0342,0.0323,0.0235,0.0246};
    std::vector<Real> u = {4.,2.,1.,0.5,0.25,0.167,0.125,0.1,
      0.0833,0.0714,0.0625};

  public:

    KowalikOsborneProblem()
      : TestLeastSquaresProblem(4,11,2,1)
    {}

    std::string name() const override
    {
      return "Kowalik and Osborne function";
    }

    Real component(unsigned int i, const Point& x) const override
    {
      return y[i] - x[0]*(std::pow(u[i],2) + u[i]*x[1])/(std::pow(u[i],2) + u[i]*x[2] + x[3]);
    }

    void gradient_comp(unsigned int i, const Point& x,
        Point& gradient) const override
    {
      gradient[0] += -2. * (std::pow(u[i],2) + u[i]*x[1])/(std::pow(u[i],2) + u[i]*x[2] + x[3]) * component(i,x);
      gradient[1] += -2. * u[i]*x[0]/(std::pow(u[i],2) + u[i]*x[2] + x[3]) * component(i,x);
      gradient[2] +=  2. * u[i]*x[0]*(std::pow(u[i],2) + u[i]*x[1])/std::pow(std::pow(u[i],2) + u[i]*x[2] + x[3],2)
        * component(i,x);
      gradient[3] +=  2. * x[0]*(std::pow(u[i],2) + u[i]*x[1])/std::pow(std::pow(u[i],2) + u[i]*x[2] + x[3],2)
        * component(i,x);
    }

    Point start_point(Real scale = 1.) const override
    {
      return {0.25*scale,0.39*scale,0.415*scale,0.39*scale};
    }

    Real goal_value(unsigned int i) const override
    {
      if (i == 0)
        return 3.075056038512582e-4;
      else
        return 1.02734e-3;
    }

    Point goal_point(unsigned int i) const override
    {
      // unknown, best found point
      return {0.192807,0.191282,0.123057,0.136062};
    }
};

#endif // DUNE_NONLINOPT_TEST_KOWALIKOSBORNE_HH
