#ifndef DUNE_NONLINOPT_TEST_JENNRICHSAMPSON_HH
#define DUNE_NONLINOPT_TEST_JENNRICHSAMPSON_HH

#include "leastsquares.hh"

/**
 * @brief Test problem 6 from Moré, Garbow and Hillstrom
 */
class JennrichSampsonProblem
: public TestLeastSquaresProblem
{
  public:

    using Real  = typename TestLeastSquaresProblem::Real;
    using Point = typename TestLeastSquaresProblem::Point;

    JennrichSampsonProblem(unsigned int m_ = 10)
      : TestLeastSquaresProblem(2,m_)
    {
      if (m_ < 2)
      {
#if HAVE_DUNE_COMMON
        DUNE_THROW(Dune::Exception,"m >= 2 required");
#else
        throw std::logic_error("m >= 2 required");
#endif // HAVE_DUNE_COMMON
      }
    }

    std::string name() const override
    {
      return "Jennrich and Sampson function";
    }

    Real component(unsigned int i, const Point& x) const override
    {
      return 2. + 2.*(i+1) - (std::exp((i+1)*x[0]) + std::exp((i+1)*x[1]));
    }

    void gradient_comp(unsigned int i, const Point& x,
        Point& gradient) const override
    {
      gradient[0] += -2. * (i+1)*std::exp((i+1)*x[0])
        * (2. + 2.*(i+1) - (std::exp((i+1)*x[0]) + std::exp((i+1)*x[1])));
      gradient[1] += -2. * (i+1)*std::exp((i+1)*x[1])
        * (2. + 2.*(i+1) - (std::exp((i+1)*x[0]) + std::exp((i+1)*x[1])));
    }

    Point start_point(Real scale = 1.) const override
    {
      return {0.3*scale,0.4*scale};
    }

    Real goal_value(unsigned int i) const override
    {
      return 1.243621823556147e+02;
    }

    Point goal_point(unsigned int i) const override
    {
      return {0.257825,0.257825};
    }
};

#endif // DUNE_NONLINOPT_TEST_JENNRICHSAMPSON_HH
