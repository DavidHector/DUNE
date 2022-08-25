#ifndef DUNE_NONLINOPT_TEST_ROSENBROCK_HH
#define DUNE_NONLINOPT_TEST_ROSENBROCK_HH

#include "leastsquares.hh"

/**
 * @brief Test problem 1 from Moré, Garbow and Hillstrom
 */
class RosenbrockProblem
: public TestLeastSquaresProblem
{
  public:

    using Real  = typename TestLeastSquaresProblem::Real;
    using Point = typename TestLeastSquaresProblem::Point;

  private:

    double a,b;

  public:

    RosenbrockProblem(double a_ = 1., double b_ = 10., unsigned int n_ = 2)
      : TestLeastSquaresProblem(n_,n_), a(a_), b(b_)
    {
      if (this->n % 2 != 0)
      {
#if HAVE_DUNE_COMMON
        DUNE_THROW(Dune::Exception,"n must be even");
#else
        throw std::logic_error("n must be even");
#endif // HAVE_DUNE_COMMON
      }
    }

    std::string name() const override
    {
      if (this->n == 2)
        return "Rosenbrock function";
      else
        return "Extended Rosenbrock function";
    }

    Real component(unsigned int i, const Point& x) const override
    {
      if (i % 2 == 0)
        return x[i] - a;
      else
        return b * (std::pow(x[i-1],2) - x[i]);
    }

    void gradient_comp(unsigned int i, const Point& x,
        Point& gradient) const override
    {
      if (i % 2 == 0)
        gradient[i] += -2. * (a - x[i]);
      else
      {
        const Real factor = std::pow(b,2) * (x[i] - std::pow(x[i-1],2));
        gradient[i-1] += -4. * x[i-1] * factor;
        gradient[i] += 2. * factor;
      }
    }

    Point start_point(Real scale = 1.) const override
    {
      Point start = zero();

      for (unsigned int i = 0; i < n; i += 2)
      {
        start[i] = scale * (-1.2 + i*0.1);
        start[i+1] = scale;
      }

      return start;
    }

    Real goal_value(unsigned int i) const override
    {
      return 0.;
    }

    Point goal_point(unsigned int i) const override
    {
      return Point(n,1.);
    }
};

#endif // DUNE_NONLINOPT_TEST_ROSENBROCK_HH
