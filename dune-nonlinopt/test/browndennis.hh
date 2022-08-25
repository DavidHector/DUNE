#ifndef DUNE_NONLINOPT_TEST_BROWNDENNIS_HH
#define DUNE_NONLINOPT_TEST_BROWNDENNIS_HH

#include "leastsquares.hh"

/**
 * @brief Test problem 16 from Moré, Garbow and Hillstrom
 */
class BrownDennisProblem
: public TestLeastSquaresProblem
{
  public:

    using Real  = typename TestLeastSquaresProblem::Real;
    using Point = typename TestLeastSquaresProblem::Point;

    BrownDennisProblem(unsigned int m_ = 20)
      : TestLeastSquaresProblem(4,m_,1,1)
    {
      if (m_ < 4)
      {
#if HAVE_DUNE_COMMON
        DUNE_THROW(Dune::Exception,"m must be at least 4");
#else
        throw std::logic_error("m must be at least 4");
#endif // HAVE_DUNE_COMMON
      }
    }

    std::string name() const override
    {
      return "Brown and Dennis function";
    }

    Real component(unsigned int i, const Point& x) const override
    {
      const Real t_i = (i+1)/5.;
      return std::pow(x[0] + t_i*x[1] - std::exp(t_i),2)
        + std::pow(x[2] + x[3]*std::sin(t_i) - std::cos(t_i),2);
    }

    void gradient_comp(unsigned int i, const Point& x,
        Point& gradient) const override
    {
      const Real t_i = (i+1)/5.;
      gradient[0] += 2. * 2.*(x[0] + t_i*x[1] - std::exp(t_i)) * component(i,x);
      gradient[1] += 2. * 2.*t_i*(x[0] + t_i*x[1] - std::exp(t_i)) * component(i,x);
      gradient[2] += 2. * 2.*(x[2] + x[3]*std::sin(t_i) - std::cos(t_i)) * component(i,x);
      gradient[3] += 2. * 2.*std::sin(t_i)
        * (x[2] + x[3]*std::sin(t_i) - std::cos(t_i)) * component(i,x);
    }

    Point start_point(Real scale = 1.) const override
    {
      return {25.*scale,5.*scale,-5.*scale,-scale};
    }

    Real goal_value(unsigned int i) const override
    {
      return 8.582220162635627e+04;
    }

    Point goal_point(unsigned int i) const override
    {
      // unknown, best found point
      return {-11.5944,13.2036,-0.403439,0.236779};
    }
};

#endif // DUNE_NONLINOPT_TEST_BROWNDENNIS_HH
