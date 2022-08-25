#ifndef DUNE_NONLINOPT_TEST_BIGGSEXP6_HH
#define DUNE_NONLINOPT_TEST_BIGGSEXP6_HH

#include "leastsquares.hh"

/**
 * @brief Test problem 18 from Moré, Garbow and Hillstrom
 */
class BiggsEXP6Problem
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

    BiggsEXP6Problem(unsigned int m_ = 13)
      : TestLeastSquaresProblem(6,m_,2,2)
    {
      if (this->m < 6)
      {
#if HAVE_DUNE_COMMON
        DUNE_THROW(Dune::Exception,"m must be >= 6");
#else
        throw std::logic_error("m must be >= 6");
#endif // HAVE_DUNE_COMMON
      }
    }

    std::string name() const override
    {
      return "Biggs EXP6 function";
    }

    Real component(unsigned int i, const Point& x) const override
    {
      const Real t_i = 0.1*(i+1);
      const Real y_i = std::exp(-t_i) - 5.*std::exp(-10.*t_i) + 3.*std::exp(-4.*t_i);
      return x[2]*std::exp(-t_i*x[0]) - x[3]*std::exp(-t_i*x[1]) + x[5]*std::exp(-t_i*x[4]) - y_i;
    }

    void gradient_comp(unsigned int i, const Point& x,
        Point& gradient) const override
    {
      const Real t_i = 0.1*(i+1);
      gradient[0] += -2. * t_i*x[2]*std::exp(-t_i*x[0]) * component(i,x);
      gradient[1] +=  2. * t_i*x[3]*std::exp(-t_i*x[1]) * component(i,x);
      gradient[2] +=  2. * std::exp(-t_i*x[0]) * component(i,x);
      gradient[3] += -2. * std::exp(-t_i*x[1]) * component(i,x);
      gradient[4] += -2. * t_i*x[5]*std::exp(-t_i*x[4]) * component(i,x);
      gradient[5] +=  2. * std::exp(-t_i*x[4]) * component(i,x);
    }

    Point start_point(Real scale = 1.) const override
    {
      return {scale,2.*scale,scale,scale,scale,scale};
    }

    Real goal_value(unsigned int i) const override
    {
      if (i == 0)
        return 0.;
      else
        return 5.655649925505601e-03;
    }

    Point goal_point(unsigned int i) const override
    {
      if (i == 0)
        return {1.,10.,1.,5.,4.,3.};
      else
        // unknown, best found point
        return {1.71142,17.6832,1.16314,5.18656,1.71142,1.16314};
    }
};

#endif // DUNE_NONLINOPT_TEST_BIGGSEXP6_HH
