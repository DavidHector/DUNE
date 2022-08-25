#ifndef DUNE_NONLINOPT_TEST_LEASTSQUARES_HH
#define DUNE_NONLINOPT_TEST_LEASTSQUARES_HH

#include<dune/nonlinopt/vectorclass.hh>
#include<dune/nonlinopt/nonlinopt.hh>

#include<iomanip>

/**
 * @brief Interface for the test least squares problems
 *
 * This class provides a common interface for the test
 * problems of the automatic testing suite.
 */
class TestLeastSquaresProblem
: public Dune::NonlinOpt::ProblemBase<double,Dune::NonlinOpt::VectorClass<double>>
{
  public:

    using Real  = double;
    using Point = Dune::NonlinOpt::VectorClass<Real>;

  protected:

    unsigned int n,m,v,p;
    mutable unsigned int iteration = 0;
    mutable unsigned int value_count = 0, grad_count = 0;

  public:

    /**
     * @brief Constructor
     *
     * @param n_ dimension of the problem
     * @param m_ number of least squares components
     * @param v_ number of known minimal values
     * @param p_ number of known local minimizers (positions)
     */
    TestLeastSquaresProblem(unsigned int n_, unsigned int m_,
        unsigned int v_ = 1, unsigned int p_ = 1)
      : n(n_), m(m_), v(v_), p(p_)
    {}

    /**
     * @brief Name of test problem
     *
     * @return string containing the name
     */
    virtual std::string name() const = 0;

    Point zero() const override
    {
      return Point(n,0.);
    }

    std::size_t dim() const override
    {
      return n;
    }

    Real value(const Point& x, bool subsequent = false) const override
    {
      if (x.size() != n)
      {
#if HAVE_DUNE_COMMON
        DUNE_THROW(Dune::Exception,"Point has wrong dimension");
#else
        throw std::logic_error("Point has wrong dimension");
#endif // HAVE_DUNE_COMMON
      }

      Real output = 0.;
      for (unsigned int i = 0; i < m; i++)
        output += std::pow(component(i,x),2);

      value_count++;
      return output;
    }

    Real dir_deriv(const Point& x, const Point& direction,
        Point& gradient, bool subsequent = false) const override
    {
      this->gradient(x,gradient);
      return gradient * direction;
    }

    void gradient(const Point& x, Point& gradient, const Point* const direction = nullptr,
        Real derivative = 0.) const override
    {
      if (direction)
        return;

      if (x.size() != n)
      {
#if HAVE_DUNE_COMMON
        DUNE_THROW(Dune::Exception,"Point has wrong dimension");
#else
        throw std::logic_error("Point has wrong dimension");
#endif // HAVE_DUNE_COMMON
      }

      gradient = zero();

      for (unsigned int i = 0; i < m; i++)
        gradient_comp(i,x,gradient);

      grad_count++;
    }

    std::pair<Real,Real> value_and_deriv(const Point& x, const Point& direction,
        Point& gradient) const override
    {
      this->gradient(x,gradient);
      return {value(x), gradient * direction};
    }

    void hook(std::size_t iter, const Point& x, Real value,
        const Point& gradient, bool extrapolation = false) const override
    {
      iteration = iter;
    }

    /**
     * @brief Contribution of i-th least squares component to value
     *
     * @param i index of component to evaluate
     * @param x vector of function variables
     *
     * @return function value
     */
    virtual Real component(unsigned int i, const Point& x) const = 0;
    
    /**
     * @brief Contribution of i-th least squares component to gradient
     *
     * @param i       index of component to evaluate
     * @param x       vector of function variables
     * @param[in,out] gradient, current contribution will be added
     */
    virtual void gradient_comp(unsigned int i, const Point& x, Point& gradient) const = 0;

    /**
     * @brief Initial vector of function variables
     *
     * @param scale scales start point to generate family of test problems
     */
    virtual Point start_point(Real scale = 1.) const = 0;

    /**
     * @brief Location of i-th known minimum
     *
     * @param i index to select among several minima
     *
     * @return vector containing location of minimum
     */
    virtual Point goal_point(unsigned int i) const = 0;

    /**
     * @brief Value of i-th known minimum
     *
     * @param i index to select among several minima
     *
     * @return i-th local minimum value
     */
    virtual Real goal_value(unsigned int i) const = 0;

    /**
     * @brief Distance of parameters from closest known minimizer
     *
     * Computes the Euclidean distance to the set of known
     * minimizers, i.e., the distance to the optimal location if
     * there is only one, or the closest minimizer if there are several.
     *
     * @param point vector of function parameters
     *
     * @return value containing Euclidean distance
     */
    Real error(const Point& point) const
    {
      Real output = std::numeric_limits<Real>::max();
      for (unsigned int i = 0; i < p; i++)
      {
        Point goal = goal_point(i);
        Real error_i = 0;
        for (unsigned int j = 0; j < n; j++)
          error_i += std::pow(point[j]-goal[j],2);
        if (error_i < output)
          output = error_i;
      }
      return std::sqrt(output);
    }

    /**
     * @brief Distance of function value from clostest known minimum
     *
     * Computes the absolute value of the difference between the
     * current function value and the optimal one, respectively the
     * smallest such difference if there are several local minima with
     * differing function values.
     *
     * @param val current function value
     *
     * @return distance to closest known local minimum
     */
    Real residual(Real val) const
    {
      Real output = std::numeric_limits<Real>::max();
      Real minAbs = std::numeric_limits<Real>::max();
      for (unsigned int i = 0; i < v; i++)
      {
        const Real minVal = goal_value(i);
        if (std::abs(val - minVal) < minAbs)
        {
          minAbs = std::abs(val - minVal);
          output = val - minVal;
        }
      }
      return output;
    }

    /**
     * @brief Convergence test based on known values
     *
     * This function test whether an optimization run that claims to
     * have converged actually reached one of the known minima. This is
     * deemed to be the case if the method has reached one of the
     * function values at the known local minima with a relative tolerance
     * of 1e-5, or an absolute tolerance of 1e-5 if the minimum happens
     * to be zero.
     *
     * @param val current function value
     *
     * @return whether known minimum value was reached
     */
    bool converged(Real val) const
    {
      Real minAbs = std::numeric_limits<Real>::max();
      Real min = 0.;
      for (unsigned int i = 0; i < v; i++)
      {
        const Real minVal = goal_value(i);
        if (std::abs(val - minVal) < minAbs)
        {
          minAbs = std::abs(val - minVal);
          min = minVal;
        }
      }
      if (min == 0.)
        return (minAbs < 1e-5);
      else
        return (minAbs/min < 1e-5);
    }

    /**
     * @brief Report information about the optimization run
     *
     * This method returns the number of iterations that were taken,
     * the number of function and gradient computations, their sum
     * and weighted sum, whether the method has converged to a known
     * local minimum, the Euclidean norm of the gradient after the
     * final iteration, the resulting residual, and the resulting error.
     *
     * @param point vector of function parameters
     *
     * @return tuple containing the aforementioned values
     *
     * @see residual
     * @see error
     */
    std::tuple<unsigned int, unsigned int, unsigned int, unsigned int,
      unsigned int, bool, Real, Real, Real>
      report(const Point& point, Real val) const
    {
      Real err = error(point);
      Real res = residual(val);
      Point gradient = zero();
      this->gradient(point,gradient);

      return {iteration, value_count, grad_count, value_count + grad_count,
        value_count + 3*grad_count, converged(val), std::sqrt(gradient * gradient),
        res, err};
    }

    /**
     * @brief Check known local minima for consistency
     *
     * This method calculates the gradients at each of the local minima
     * specified in the test problem definition. The gradient should be
     * zero at these points, and non-zero gradients imply some form of
     * bug in the implementation of the test problem.
     */
    void test_values() const
    {
      for (unsigned int i = 0; i < p; i++)
      {
        Point gradient = zero();
        this->gradient(goal_point(i),gradient);
        std::cout << "goal value " << i << ": " << goal_value(i)
          << " " << value(goal_point(i)) << " gradient:";
        for (unsigned int j = 0; j < n; j++)
          std::cout << " " << gradient[j];
        std::cout << std::endl;
      }
    }
};

/**
 * @brief Function representing an optimization test run
 *
 * This version of the test run uses the default configuration.
 *
 * @param problem  instance of test problem to solve
 * @param scale    scale to apply to starting position vector
 * @param max_iter maximum number of iterations to perform
 *
 * @return information about run, see report method of problem class
 */
std::tuple<unsigned int, unsigned int, unsigned int, unsigned int,
  unsigned int, bool, double, double, double>
  solve(const TestLeastSquaresProblem& problem, double scale, unsigned int max_iter = 9999);

#if HAVE_DUNE_COMMON
/**
 * @brief Function representing an optimization test run
 *
 * This version of the test run can be configured through
 * a Dune::ParameterTree object.
 *
 * @param problem  instance of test problem to solve
 * @param scale    scale to apply to starting position vector
 * @param config   solver configuration object
 *
 * @return information about run, see report method of problem class
 */
std::tuple<unsigned int, unsigned int, unsigned int, unsigned int,
  unsigned int, bool, double, double, double>
  solve(const TestLeastSquaresProblem& problem, double scale,
      const Dune::ParameterTree& config);
#endif // HAVE_DUNE_COMMON

#endif // DUNE_NONLINOPT_TEST_LEASTSQUARES_HH
