#ifndef DUNE_NONLINOPT_TEST_BOX3D_HH
#define DUNE_NONLINOPT_TEST_BOX3D_HH

#include "leastsquares.hh"

/**
 * @brief Test problem 12 from Moré, Garbow and Hillstrom
 */
class Box3DProblem
: public TestLeastSquaresProblem
{
  public:

    using Real  = typename TestLeastSquaresProblem::Real;
    using Point = typename TestLeastSquaresProblem::Point;

    Box3DProblem(unsigned int m_ = 100)
      : TestLeastSquaresProblem(3,m_,1,2)
    {}

    std::string name() const override
    {
      return "Box three-dimensional function";
    }

    Real component(unsigned int i, const Point& x) const override
    {
      const Real t_i = 0.1*(i+1);
      return std::exp(-t_i*x[0]) - std::exp(-t_i*x[1])
        - x[2]*(std::exp(-t_i) - std::exp(-10.*t_i));
    }

    void gradient_comp(unsigned int i, const Point& x,
        Point& gradient) const override
    {
      const Real t_i = 0.1*(i+1);
      gradient[0] += -2. * t_i*std::exp(-t_i*x[0]) * component(i,x);
      gradient[1] +=  2. * t_i*std::exp(-t_i*x[1]) * component(i,x);
      gradient[2] +=  2. * (std::exp(-10.*t_i) - std::exp(-t_i)) * component(i,x);
    }

    Point start_point(Real scale = 1.) const override
    {
      return {0.,10.*scale,20.*scale};
    }

    Real goal_value(unsigned int i) const override
    {
      return 0.;
    }

    Point goal_point(unsigned int i) const override
    {
      // several, (1,10,1), (10,1,-1), (x,x,0)
      if (i == 0)
        return {1.,10.,1.};
      else
        return {10.,1.,-1.};
    }
};

#endif // DUNE_NONLINOPT_TEST_BOX3D_HH
