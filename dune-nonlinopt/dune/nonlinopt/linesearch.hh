#ifndef DUNE_NONLINOPT_LINESEARCH_HH
#define DUNE_NONLINOPT_LINESEARCH_HH

#include<iostream>
#include<iomanip>

#if HAVE_DUNE_COMMON
#include<dune/common/parametertree.hh>
#endif // HAVE_DUNE_COMMON

#include<dune/nonlinopt/linesearchcriterion.hh>

namespace Dune {
  namespace NonlinOpt {

    //! @brief General exception for line search failures
#if HAVE_DUNE_COMMON
    struct LinesearchFailed    : public Dune::Exception {};
#else
    struct LinesearchFailed    : public std::exception {};
#endif // HAVE_DUNE_COMMON

    //! @brief Directional derivative is non-negative
    struct NoDescentDirection  : public LinesearchFailed {};
    //! @brief Rapid expansion phase failed
    struct SlopeAlwaysNegative : public LinesearchFailed {};
    //! @brief Bracketing during contraction phase failed
    struct BracketingFailed    : public LinesearchFailed {};

    template<typename,typename> class ProblemBase;

    /**
     * @brief Abstract base class for line search methods
     *
     * Derive from this class to provide a custom line search.
     */
    template<typename Real, typename Point>
      class LinesearchBase
      {
        protected:

          //! @brief Line search termination criterion
          std::unique_ptr<LinesearchCriterionBase<Real,Point>> criterion;

        public:

          /**
           * @brief Default constructor
           *
           * Sets relaxed Wolfe conditions as line search termination
           * criterion. Other criteria can be selected through the
           * corresponding setter method.
           *
           * @see set_linesearchcriterion
           */
          LinesearchBase()
          {
            set_linesearchcriterion<RelaxedWolfeConditions>();
          }

#if HAVE_DUNE_COMMON
          /**
           * @brief Constructor based on Dune::ParameterTree
           *
           * Chooses line search termination criterion based on
           * configuration, default is relaxed Wolfe conditions.
           *
           * @param config ParameterTree object with configuration
           */
          LinesearchBase(const Dune::ParameterTree& config)
          {
            const std::string crit_name = config.get<std::string>("linesearch.criterion","nonmonotone_relaxed_wolfe");
            if (crit_name == "nonmonotone_relaxed_wolfe")
              criterion = std::make_unique<NonmonotoneRelaxedWolfeConditions<Real,Point>>(config);
            else if (crit_name == "relaxed_wolfe")
              criterion = std::make_unique<RelaxedWolfeConditions<Real,Point>>(config);
            else if (crit_name == "approx_wolfe")
              criterion = std::make_unique<ApproximateWolfeConditions<Real,Point>>(config);
            else if (crit_name == "wolfe")
              criterion = std::make_unique<WolfeConditions<Real,Point>>(config);
            else if (crit_name == "strong_wolfe")
              criterion = std::make_unique<StrongWolfeConditions<Real,Point>>(config);
            else if (crit_name == "nonmonotone_armijo")
              criterion = std::make_unique<NonmonotoneArmijoRule<Real,Point>>(config);
            else if (crit_name == "armijo")
              criterion = std::make_unique<ArmijoRule<Real,Point>>(config);
            else if (crit_name == "goldstein")
              criterion = std::make_unique<GoldsteinConditions<Real,Point>>(config);
            else
              DUNE_THROW(Dune::Exception,"Linesearch criterion not known, use set_linesearchcriterion() method");
          }
#endif // HAVE_DUNE_COMMON

          /**
           * @brief Provide custom line search termination criterion
           *
           * @param args arguments passed to criterion constructor
           */
          template<template<typename,typename> class Criterion, typename... Args>
            void set_linesearchcriterion(Args&&... args)
            {
              criterion = std::make_unique<Criterion<Real,Point>>(std::forward<Args>(args)...);
            }

          /**
           * @brief Unique identifier for the line search
           *
           * @return string containing the name
           */
          virtual std::string name() const = 0;

          /**
           * @brief Unique identifier for the line search termination criterion
           *
           * @return string containing the name
           */
          std::string criterion_name()
          {
            return criterion->name();
          }

          /**
           * @brief Reset all information, for starting a new run
           *
           * Line search methods that don't have internal state
           * can simply implement a function that forwards the call
           * to the line search criterion.
           */
          virtual void hard_reset() = 0;

          /**
           * @brief Perform line search in given direction
           *
           * Applies a line search strategy along the given direction
           * until a point is found that is acceptable according to the
           * line search termination criterion. Updates its arguments
           * with the accepted step width, the resulting point, its
           * associated function value and gradient, and the directional
           * derivative at that point along the line.
           *
           * @param         problem     problem definition, for function values and gradients
           * @param         direction   current search direction
           * @param         derivative  directional derivative at start of line
           * @param[in,out] point       current point, becomes updated point
           * @param[in,out] value       current value, becomes updated value
           * @param[in,out] gradient    current gradient, becomes updated gradient
           * @param[in,out] alpha       initial guess, becomes accepted step width
           * @param[out]    deriv_alpha directional derivative at accepted step width
           */
          virtual void apply(const ProblemBase<Real,Point>& problem, const Point& direction,
              Real derivative, Point& point, Real& value, Point& gradient,
              Real& alpha, Real& deriv_alpha) const = 0;

          virtual ~LinesearchBase(){};
      };



    /**
     * @brief Simple backtracking line search
     *
     * This line search starts at the proposal step width and
     * repeatedly multiplies it with a small constant until the
     * termination criterion is fulfilled. By construction, this
     * simple line search is unable to increase the step width
     * beyond the initial proposal, which means this starting
     * value has to be chosen large enough.
     */
    template<typename Real, typename Point>
      class BacktrackingLinesearch
      : public LinesearchBase<Real,Point>
      {
        const Real        reduction;
        const std::size_t linesearch_tries;

        mutable Point test_point, test_grad;

        public:

        /**
         * @brief Constructor
         *
         * @param reduction_        factor used in backtracking
         * @param linesearch_tries_ maximum number of tested step widths
         */
        BacktrackingLinesearch(Real reduction_ = 0.9,
            std::size_t linesearch_tries_ = 100)
          : reduction(reduction_), linesearch_tries(linesearch_tries_)
        {
          if (reduction <= 0. || reduction >= 1.)
            throw std::invalid_argument("reduction has to be in (0,1)");
          if (linesearch_tries == 0)
            throw std::invalid_argument("number of tries for line search must be positive");
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts reduction factor and number of tries from configuration.
         *
         * @param config ParameterTree configuration object
         */
        BacktrackingLinesearch(const Dune::ParameterTree& config)
          : LinesearchBase<Real,Point>(config),
          reduction       (config.get<Real>       ("linesearch.backtracking_factor",0.9)),
          linesearch_tries(config.get<std::size_t>("linesearch.linesearch_tries",   100))
        {
          if (reduction <= 0. || reduction >= 1.)
            DUNE_THROW(Dune::Exception,"reduction has to be in (0,1)");
          if (linesearch_tries == 0)
            DUNE_THROW(Dune::Exception,"number of tries for line search must be positive");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "backtracking";
        }

        void hard_reset() override
        {
          this->criterion->hard_reset();
        }

        void apply(const ProblemBase<Real,Point>& problem, const Point& direction, Real derivative,
            Point& point, Real& value, Point& gradient, Real& alpha, Real& deriv_alpha) const override
        {
          std::size_t iter = 1;

          while(true)
          {
            test_point = point;
            test_point.axpy(direction,alpha);
            const Real test_value = problem.value(test_point);

            if (this->criterion->check_value(value,derivative,alpha,test_value))
            {
              const Real test_deriv = problem.dir_deriv(test_point,direction,test_grad,true);

              if (this->criterion->check_derivative(value,derivative,alpha,test_deriv))
              {
                problem.gradient(test_point,test_grad,&direction,test_deriv);

                std::swap(point,test_point);
                value = test_value;
                std::swap(gradient,test_grad);
                deriv_alpha = test_deriv;
                return;
              }
            }
            else
              alpha *= reduction;

            iter++;
            if (iter > linesearch_tries)
            {
#if HAVE_DUNE_COMMON
              DUNE_THROW(LinesearchFailed,"backtracking linesearch failed");
#else
              throw LinesearchFailed();
#endif // HAVE_DUNE_COMMON
            }
          }
        }

      };



    /**
     * @brief Simple line search based on cubic interpolation
     *
     * This line search performs an inital quadratic interpolation
     * of the function along the line, and then uses the last two
     * considered step widths with resulting function values and
     * directional derivatives for a cubic Hermite interpolation to
     * obtain a new candidate step width.
     */
    template<typename Real, typename Point>
      class CubicLinesearch
      : public LinesearchBase<Real,Point>
      {
        const std::size_t linesearch_tries;

        mutable Point test_point, test_grad;

        public:

        /**
         * @brief Constructor
         *
         * @param linesearch_tries_ maximum number of tested step widths
         */
        CubicLinesearch(std::size_t linesearch_tries_ = 100)
          : linesearch_tries(linesearch_tries_)
        {
          if (linesearch_tries == 0)
            throw std::invalid_argument("number of tries for line search must be positive");
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts number of tries from configuration.
         *
         * @param config ParameterTree configuration object
         */
        CubicLinesearch(const Dune::ParameterTree& config)
          : LinesearchBase<Real,Point>(config),
          linesearch_tries(config.get<std::size_t>("linesearch.linesearch_tries",100))
        {
          if (linesearch_tries == 0)
            DUNE_THROW(Dune::Exception,"number of tries for line search must be positive");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "cubic";
        }

        void hard_reset() override
        {
          this->criterion->hard_reset();
        }

        void apply(const ProblemBase<Real,Point>& problem, const Point& direction, Real derivative,
            Point& point, Real& value, Point& gradient, Real& alpha, Real& deriv_alpha) const override
        {
          std::size_t iter = 1;

          test_point = point;
          test_point.axpy(direction,alpha);
          Real test_value = problem.value(test_point);

          if (this->criterion->check_value(value,derivative,alpha,test_value))
          {
            const Real test_deriv = problem.dir_deriv(test_point,direction,test_grad,true);

            if (this->criterion->check_derivative(value,derivative,alpha,test_deriv))
            {
              problem.gradient(test_point,test_grad,&direction,test_deriv);

              std::swap(point,test_point);
              value = test_value;
              std::swap(gradient,test_grad);
              deriv_alpha = test_deriv;
              return;
            }
          }

          // quadratic interpolation
          Real prev_alpha = alpha;
          Real prev_test_value = test_value;
          alpha = - derivative * alpha * alpha / (2. * (test_value - value - derivative * alpha));

          while(true)
          {
            test_point = point;
            test_point.axpy(direction,alpha);
            test_value = problem.value(test_point);

            if (this->criterion->check_value(value,derivative,alpha,test_value))
            {
              const Real test_deriv = problem.dir_deriv(test_point,direction,test_grad,true);

              if (this->criterion->check_derivative(value,derivative,alpha,test_deriv))
              {
                problem.gradient(test_point,test_grad,&direction,test_deriv);

                std::swap(point,test_point);
                value = test_value;
                std::swap(gradient,test_grad);
                deriv_alpha = test_deriv;
                return;
              }
            }

            // cubic interpolation
            const Real factor = 1./(alpha*alpha * prev_alpha*prev_alpha * (prev_alpha - alpha));
            const Real y0 = test_value      - value - derivative * alpha;
            const Real y1 = prev_test_value - value - derivative * prev_alpha;
            const Real a = factor * (  alpha*alpha       * y1 - prev_alpha*prev_alpha            * y0);
            const Real b = factor * (- alpha*alpha*alpha * y1 + prev_alpha*prev_alpha*prev_alpha * y0);
            prev_alpha = alpha;
            prev_test_value = test_value;
            if (3.*a * derivative < b*b)
              alpha = (- b + std::sqrt(b*b - 3.*a * derivative)) / (3.*a);
            else
              alpha = -b / (3.*a);

            iter++;
            if (iter > linesearch_tries)
            {
#if HAVE_DUNE_COMMON
              DUNE_THROW(LinesearchFailed,"cubic line search failed");
#else
              throw LinesearchFailed();
#endif // HAVE_DUNE_COMMON
            }
          }
        }
      };



    /**
     * @brief Simple quadratic line search
     *
     * This line search uses repeated quadratic interpolation
     * to search for an acceptable step width. In contrast to
     * the cubic line search, information from previous tries
     * is not reused and two function values need to be obtained
     * instead. Therefore, a larger number of function evaluations
     * may be needed.
     *
     * @see CubicLinesearch
     */
    template<typename Real, typename Point>
      class QuadraticLinesearch
      : public LinesearchBase<Real,Point>
      {
        const std::size_t linesearch_tries;

        mutable Point test_point, test_grad;

        public:

        /**
         * @brief Constructor
         *
         * @param linesearch_tries_ maximum number of tested step widths
         */
        QuadraticLinesearch(std::size_t linesearch_tries_ = 100)
          : linesearch_tries(linesearch_tries_)
        {
          if (linesearch_tries == 0)
            throw std::invalid_argument("number of tries for line search must be positive");
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts number of tries from configuration.
         *
         * @param config ParameterTree configuration object
         */
        QuadraticLinesearch(const Dune::ParameterTree& config)
          : LinesearchBase<Real,Point>(config),
          linesearch_tries(config.get<std::size_t>("linesearch.linesearch_tries",100))
        {
          if (linesearch_tries == 0)
            DUNE_THROW(Dune::Exception,"number of tries for line search must be positive");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "quadratic";
        }

        void hard_reset() override
        {
          this->criterion->hard_reset();
        }

        void apply(const ProblemBase<Real,Point>& problem, const Point& direction, Real derivative,
            Point& point, Real& value, Point& gradient, Real& alpha, Real& deriv_alpha) const override
        {
          std::size_t iter = 1;

          while(true)
          {
            // full step width
            test_point = point;
            test_point.axpy(direction,alpha);
            Real test_value = problem.value(test_point);

            if (this->criterion->check_value(value,derivative,alpha,test_value))
            {
              const Real test_deriv = problem.dir_deriv(test_point,direction,test_grad,true);

              if (this->criterion->check_derivative(value,derivative,alpha,test_deriv))
              {
                problem.gradient(test_point,test_grad,&direction,test_deriv);

                std::swap(point,test_point);
                value = test_value;
                std::swap(gradient,test_grad);
                deriv_alpha = test_deriv;
                return;
              }
            }

            // half step width
            test_point.axpy(direction,-0.5*alpha);
            Real half_test_value = problem.value(test_point);

            if (this->criterion->check_value(value,derivative,alpha,test_value))
            {
              const Real test_deriv = problem.dir_deriv(test_point,direction,test_grad,true);

              if (this->criterion->check_derivative(value,derivative,alpha,test_deriv))
              {
                problem.gradient(test_point,test_grad,&direction,test_deriv);

                std::swap(point,test_point);
                value = test_value;
                std::swap(gradient,test_grad);
                deriv_alpha = test_deriv;
                return;
              }
            }

            // quadratic interpolation
            Real c = value;
            Real b = -(test_value - 4*half_test_value + 3*value)/alpha;
            Real a = 2*(test_value - 2*half_test_value + value)/(alpha*alpha);

            if (a > 0.)
              alpha = std::max(0.1*alpha, std::min(10.*alpha, -0.5 * b/a));
            else
              alpha = std::max(0.1*alpha, std::min(10.*alpha, -0.5 * (b + std::sqrt(b*b - 4.*a*c))/a));

            iter++;
            if (iter > linesearch_tries)
            {
#if HAVE_DUNE_COMMON
              DUNE_THROW(LinesearchFailed,"quadratic line search failed");
#else
              throw LinesearchFailed();
#endif // HAVE_DUNE_COMMON
            }
          }
        }

      };



    /**
     * @brief Variant of simple quadratic line search
     *
     * This line search operates in the same way as the
     * quadratic line search, but the quadratic interpolation
     * is not based on the old function value and two function
     * value samples, but uses just one function value sample
     * and the old function value and directional derivative
     * instead. This makes the individual tries significantly
     * cheaper, but also less robust.
     */
    template<typename Real, typename Point>
      class DerivativeQuadraticLinesearch
      : public LinesearchBase<Real,Point>
      {
        const std::size_t linesearch_tries;

        mutable Point test_point, test_grad;

        public:

        /**
         * @brief Constructor
         *
         * @param linesearch_tries_ maximum number of tested step widths
         */
        DerivativeQuadraticLinesearch(std::size_t linesearch_tries_ = 100)
          : linesearch_tries(linesearch_tries_)
        {
          if (linesearch_tries == 0)
            throw std::invalid_argument("number of tries for line search must be positive");
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts number of tries from configuration.
         *
         * @param config ParameterTree configuration object
         */
        DerivativeQuadraticLinesearch(const Dune::ParameterTree& config)
          : LinesearchBase<Real,Point>(config),
          linesearch_tries(config.get<std::size_t>("linesearch.linesearch_tries",100))
        {
          if (linesearch_tries == 0)
            DUNE_THROW(Dune::Exception,"number of tries for line search must be positive");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "derivative-quadratic";
        }

        void hard_reset() override
        {
          this->criterion->hard_reset();
        }

        void apply(const ProblemBase<Real,Point>& problem, const Point& direction, Real derivative,
            Point& point, Real& value, Point& gradient, Real& alpha, Real& deriv_alpha) const override
        {
          std::size_t iter = 1;

          while(true)
          {
            test_point = point;
            test_point.axpy(direction,alpha);
            Real test_value = problem.value(test_point);

            if (this->criterion->check_value(value,derivative,alpha,test_value))
            {
              const Real test_deriv = problem.dir_deriv(test_point,direction,test_grad,true);

              if (this->criterion->check_derivative(value,derivative,alpha,test_deriv))
              {
                problem.gradient(test_point,test_grad,&direction);

                std::swap(point,test_point);
                value = test_value;
                std::swap(gradient,test_grad);
                deriv_alpha = test_deriv;
                return;
              }
            }

            Real a = (test_value - value - derivative * alpha) / (alpha * alpha);
            alpha = std::max(0.1*alpha, std::min(10.*alpha, - 0.5 * derivative / a));

            iter++;
            if (iter > linesearch_tries)
            {
#if HAVE_DUNE_COMMON
              DUNE_THROW(LinesearchFailed,"deriv. quadratic line search failed");
#else
              throw LinesearchFailed();
#endif // HAVE_DUNE_COMMON
            }
          }
        }

      };



    /**
     * @brief Simple secant line search
     *
     * This is a third variant of the quadratic line search,
     * performing a linear interpolation of the directional
     * derivative instead of a quadratic interpolation of the
     * function values. Consequently, one gradient evaluation
     * is needed per try, instead of one resp. two function
     * values.
     */
    template<typename Real, typename Point>
      class SecantLinesearch
      : public LinesearchBase<Real,Point>
      {
        const std::size_t linesearch_tries;

        mutable Point test_point, test_grad;

        public:

        /**
         * @brief Constructor
         *
         * @param linesearch_tries_ maximum number of tested step widths
         */
        SecantLinesearch(std::size_t linesearch_tries_ = 100)
          : linesearch_tries(linesearch_tries_)
        {
          if (linesearch_tries == 0)
            throw std::invalid_argument("number of tries for line search must be positive");
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts number of tries from configuration.
         *
         * @param config ParameterTree configuration object
         */
        SecantLinesearch(const Dune::ParameterTree& config)
          : LinesearchBase<Real,Point>(config),
          linesearch_tries(config.get<std::size_t>("linesearch.linesearch_tries",100))
        {
          if (linesearch_tries == 0)
            DUNE_THROW(Dune::Exception,"number of tries for line search must be positive");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "secant";
        }

        void hard_reset() override
        {
          this->criterion->hard_reset();
        }

        void apply(const ProblemBase<Real,Point>& problem, const Point& direction, Real derivative,
            Point& point, Real& value, Point& gradient, Real& alpha, Real& deriv_alpha) const override
        {
          std::size_t iter = 1;

          while(true)
          {
            test_point = point;
            test_point.axpy(direction,alpha);
            const Real test_deriv = problem.dir_deriv(test_point,direction,test_grad);

            if (this->criterion->check_derivative(value,derivative,alpha,test_deriv))
            {
              const Real test_value = problem.value(test_point,true);

              if (this->criterion->check(value,derivative,alpha,test_value,test_deriv))
              {
                problem.gradient(test_point,test_grad,&direction,test_deriv);

                std::swap(point,test_point);
                value = test_value;
                std::swap(gradient,test_grad);
                deriv_alpha = test_deriv;
                return;
              }
            }

            alpha = std::max(0.1*alpha, std::min(10.*alpha,
                  - derivative / (test_deriv - derivative) * alpha));
            iter++;
            if (iter > linesearch_tries)
            {
#if HAVE_DUNE_COMMON
              DUNE_THROW(LinesearchFailed,"gradient interp. line search failed");
#else
              throw LinesearchFailed();
#endif // HAVE_DUNE_COMMON
            }
          }

        }

      };



    /**
     * @brief Hager & Zhang's CG_DESCENT line search routine
     *
     * This is a sophisticated line search that is very robust
     * and at the same time has low average cost in terms of
     * function and gradient evaluations.
     *
     * The line search starts with a quadratic interpolation step,
     * followed by rapid expansion of the resulting search
     * interval until a step width resulting in positive slope is
     * found. Then the search interval is contracted again, while
     * guaranteeing that the interval always contains at least one
     * acceptable point.
     *
     * The original routine uses a special secant^2 step to
     * contract the interval, consisting in a pair of secant
     * interpolation and extrapolation steps. This version uses
     * cubic Hermite interpolation on the interval instead,
     * thereby making use of function values if they are available.
     * The original secant^2 step is used when function values are
     * unavailable, and the cubic interpolation can be deactivated
     * through a corresponding configuration option.
     */
    template<typename Real, typename Point>
      class HagerZhangLinesearch
      : public LinesearchBase<Real,Point>
      {
        std::size_t verbosity;
        std::size_t linesearch_tries;
        bool        combined_eval;
        bool        quad_interp;
        bool        use_cubic;
        Real        quad_threshold;
        Real        cubic_threshold;
        Real        min_expand;
        Real        max_expand;
        Real        theta;
        Real        gamma;
        Real        epsilon;

        const Real invalid = std::numeric_limits<Real>::max();

        mutable Point test_point, test_grad;
        mutable Real test_value, test_deriv;

        public:

        /**
         * @brief Constructor
         *
         * @param use_cubic_        use cubic Hermite interpolation instead of secant
         * @param combined_eval_    always evaluate values and gradients together
         * @param verbosity_        controls amount of information that is printed
         * @param linesearch_tries_ maximum iterations during expansion / contraction
         * @param quad_interp_      start with quadratic interpolation step
         * @param quad_threshold_   threshold for testing auxiliary point of interpolation
         * @param min_expand_       minimum increase factor during expansion
         * @param max_expand_       maximum increase factor during expansion
         * @param theta_            interval cutting point, defaults to center (bisection)
         * @param gamma_            threshold for performing an additional bisection step
         * @param epsilon_          estimate for relative function value error
         * @param cubic_threshold_  threshold to use vertex point if there are no roots
         */
        HagerZhangLinesearch(bool use_cubic_ = true, bool combined_eval_ = false,
            std::size_t verbosity_ = 0, std::size_t linesearch_tries_ = 100, bool quad_interp_ = true,
            Real quad_threshold_ = 0.1, Real min_expand_ = 10., Real max_expand_ = 1e6,
            Real theta_ = 0.5, Real gamma_ = 2./3., Real epsilon_ = 1e-6, Real cubic_threshold_ = 0.1)
          :
            verbosity(verbosity_), linesearch_tries(linesearch_tries_), combined_eval(combined_eval_),
            quad_interp(quad_interp_), use_cubic(use_cubic_), quad_threshold(quad_threshold_),
            cubic_threshold(cubic_threshold_), min_expand(min_expand_), max_expand(max_expand_),
            theta(theta_), gamma(gamma_), epsilon(epsilon_)
        {
          if (quad_threshold <= 0.)
            throw std::invalid_argument("quad_threshold must be positive");
          if (cubic_threshold <= 0.)
            throw std::invalid_argument("cubic_threshold must be positive");
          if (min_expand <= 0.)
            throw std::invalid_argument("min_expand must be positive");
          if (max_expand <= min_expand)
            throw std::invalid_argument("max_expand can't be smaller than min_expand");
          if (theta <= 0. || theta >= 1.)
            throw std::invalid_argument("theta must be in (0,1)");
          if (gamma <= 0. || gamma >= 1.)
            throw std::invalid_argument("gamma must be in (0,1)");
          if (epsilon < 0.)
            throw std::invalid_argument("epsilon must be non-negative");
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts parameters from configuration object.
         *
         * @param config ParameterTree configuration object
         */
        HagerZhangLinesearch(const Dune::ParameterTree& config)
          : LinesearchBase<Real,Point>(config),
          verbosity           (config.get<std::size_t>("linesearch.verbosity",         0)),
          linesearch_tries    (config.get<std::size_t>("linesearch.linesearch_tries",  100)),
          combined_eval       (config.get<bool>       ("linesearch.combined_eval",     false)),
          quad_interp         (config.get<bool>       ("linesearch.quad_interpolation",true)),
          use_cubic           (config.get<bool>       ("linesearch.use_cubic",         true)),
          quad_threshold      (config.get<Real>       ("linesearch.quad_threshold",    0.1)),
          cubic_threshold     (config.get<Real>       ("linesearch.cubic_threshold",   0.1)),
          min_expand          (config.get<Real>       ("linesearch.min_expand",        10.)),
          max_expand          (config.get<Real>       ("linesearch.max_expand",        1e6)),
          theta               (config.get<Real>       ("linesearch.hz_theta",          0.5)),
          gamma               (config.get<Real>       ("linesearch.hz_gamma",          2./3.)),
          epsilon             (config.get<Real>       ("linesearch.hz_epsilon",        1e-6))
        {
          if (quad_threshold <= 0.)
            DUNE_THROW(Dune::Exception,"quad_threshold must be positive");
          if (cubic_threshold <= 0.)
            DUNE_THROW(Dune::Exception,"cubic_threshold must be positive");
          if (min_expand <= 0.)
            DUNE_THROW(Dune::Exception,"min_expand must be positive");
          if (max_expand <= min_expand)
            DUNE_THROW(Dune::Exception,"max_expand can't be smaller than min_expand");
          if (theta <= 0. || theta >= 1.)
            DUNE_THROW(Dune::Exception,"theta must be in (0,1)");
          if (gamma <= 0. || gamma >= 1.)
            DUNE_THROW(Dune::Exception,"gamma must be in (0,1)");
          if (epsilon < 0.)
            DUNE_THROW(Dune::Exception,"epsilon must be non-negative");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "Hager-Zhang";
        }

        void hard_reset() override
        {
          this->criterion->hard_reset();
        }

        void apply(const ProblemBase<Real,Point>& problem, const Point& direction, Real derivative,
            Point& point, Real& value, Point& gradient, Real& alpha, Real& deriv_alpha) const override
        {
          Real cur_min_expand = min_expand;
          Real cur_max_expand = max_expand;

          if (verbosity >= 1)
            std::cout << "linesearch: old value " << value
              << " old deriv " << derivative
              << std::endl;
          if (derivative >= 0.)
          {
#if HAVE_DUNE_COMMON
            DUNE_THROW(NoDescentDirection,"directional derivative is positive, no descent");
#else
            throw NoDescentDirection();
#endif // HAVE_DUNE_COMMON
          }

          test_value = invalid;
          Real prev_alpha = 0., prev_deriv = derivative;
          std::size_t count = 0;
          Real a = 0.;
          Real a_value = value;
          Real a_deriv = derivative;

          // quadratic interpolation
          bool reuse_quad = false;
          if (quad_interp)
          {
            if (verbosity >= 1)
              std::cout << "initial alpha: " << alpha << std::endl;
            test_point = point;
            test_point.axpy(direction,alpha);

            Real quad_value = problem.value(test_point);

            if (std::isfinite(quad_value))
            {
              const Real denom = (quad_value - value - derivative * alpha) / std::pow(alpha,2);

              if (verbosity >= 2)
                std::cout << "interp     interval: [ " << std::setw(13) << std::scientific << a
                  << " , " << std::setw(13) << std::scientific << alpha
                  << " ] , width: " << std::setw(13) << std::scientific << alpha - a
                  << " , probe: " << std::setw(13) << std::scientific << - 0.5 * derivative / denom
                  << std::endl;

              prev_alpha = alpha;
              if (denom > 0.)
                alpha = - 0.5 * derivative / denom;

              if (prev_alpha / alpha < 1. + quad_threshold && alpha / prev_alpha < 1. + quad_threshold)
              {
                alpha = prev_alpha;
                test_value = quad_value;
                reuse_quad = true;
              }
            }
            else
              alpha /= 10.;
          }

          // rapid expansion loop
          while (true)
          {
            if (reuse_quad)
            {
              test_deriv = problem.dir_deriv(test_point,direction,test_grad,true);
              reuse_quad = false;
            }
            else
            {
              test_point = point;
              test_point.axpy(direction,alpha);

              if (combined_eval)
                std::tie(test_value, test_deriv) = problem.value_and_deriv(test_point,direction,test_grad);
              else
              {
                test_deriv = problem.dir_deriv(test_point,direction,test_grad);
                test_value = invalid;
              }
            }

            if (verbosity >= 2)
              std::cout << "expand     interval: [ " << std::setw(13) << std::scientific << a
                << " , " << std::setw(13) << std::scientific << alpha
                << " ] , width: " << std::setw(13) << std::scientific << alpha - a
                << " , probe: " << std::setw(13) << std::scientific << alpha
                << " , deriv: " << std::setw(13) << std::scientific << test_deriv << std::endl;

            if (checkCriterion(problem,direction,derivative,point,value,gradient,alpha,deriv_alpha))
            {
              return;
            }

            if (!std::isfinite(test_deriv))
            {
              cur_min_expand = std::sqrt(cur_min_expand);
              cur_max_expand = std::sqrt(cur_max_expand);
              if (verbosity >= 2)
                std::cout << "min_expand: " << cur_min_expand
                  << " max_expand: " << cur_max_expand << std::endl;
              alpha = a + 0.1 * (alpha - a);
            }
            else
            {
              if (test_deriv >= 0.)
                break;

              if (test_value <= value + epsilon * std::abs(value))
              {
                a = alpha;
                a_value = test_value;
                a_deriv = test_deriv;
              }

              const Real sec = secant(prev_alpha,prev_deriv,alpha,test_deriv);

              prev_alpha = alpha;
              prev_deriv = test_deriv;

              if (quad_interp || count != 0)
                alpha = std::max(cur_min_expand*alpha,std::min(cur_max_expand*alpha,sec));
              else
                alpha = std::max(0.1*alpha,std::min(cur_max_expand*alpha,sec));
            }

            count++;
            if (count > linesearch_tries)
            {
#if HAVE_DUNE_COMMON
              DUNE_THROW(SlopeAlwaysNegative,"slope in linesearch always negative");
#else
              throw SlopeAlwaysNegative();
#endif // HAVE_DUNE_COMMON
            }
          }

          // contraction loop
          Real b = alpha;
          Real b_value = test_value;
          Real b_deriv = test_deriv;
          while (true)
          {
            const Real old_a = a;
            const Real old_b = b;
            doubleCubic(problem,direction,derivative,point,value,gradient,
                deriv_alpha,a,a_value,a_deriv,b,b_value,b_deriv);

            if (a == b)
            {
              alpha = a;
              return;
            }
            else if (b - a > gamma * (old_b - old_a))
            {
              Real c = 0.5 * (a + b);
              updateInterval(problem,direction,derivative,point,value,gradient,
                  deriv_alpha,a,a_value,a_deriv,b,b_value,b_deriv,c);

              if (a == b)
              {
                alpha = a;
                return;
              }
            }

            count++;
            if (count > linesearch_tries)
            {
#if HAVE_DUNE_COMMON
              DUNE_THROW(LinesearchFailed,"linesearch failed to converge");
#else
              throw LinesearchFailed();
#endif // HAVE_DUNE_COMMON
            }
          }
        }

        private:

        /**
         * @brief "secant" routine of CG_DESCENT
         *
         * Performs linear interplation of the derivatives,
         * and returns the root of the interpolant.
         *
         * @param a       left interval boundary
         * @param a_deriv derivative at left boundary
         * @param b       right interval boundary
         * @param b_deriv derivative at right boundary
         *
         * @return result of interpolation resp. extrapolation
         */
        Real secant(Real a, Real a_deriv, Real b, Real b_deriv) const
        {
          const Real sec = (a * b_deriv - b* a_deriv)/(b_deriv - a_deriv);

          if (verbosity >= 2)
            std::cout << " [ " << std::setw(13) << std::scientific << a
              << " , " << std::setw(13) << std::scientific << b
              << " ] , secant: " << std::setw(13) << std::scientific
              << sec << std::endl;

          return std::isfinite(sec) ? sec : invalid;
        }

        /**
         * @brief Cubic interpolation of function, quadratic interpolation of derivative
         *
         * Performs cubic Hermite interpolation of the function values,
         * which produces a quadratic function interpolating the derivatives.
         * Uses a secant step as fallback if the function values should be
         * unavailable.
         *
         * The return value is the root that is clostest to the midpoint of
         * the interval for interpolation, and the root that is clostest to
         * the relevant interval boundary for extrapolation. The position
         * of the vertex of the corresponding parabola is returned instead
         * if the discriminant should be negative and the function value
         * at that vertex is below a configurable threshold.
         *
         * @param a       left interval boundary
         * @param a_value value at left boundary
         * @param a_deriv derivative at left boundary
         * @param b       right interval boundary
         * @param b_value value at right boundary
         * @param b_deriv derivative at right boundary
         *
         * @return result of interpolation resp. extrapolation
         */
        Real cubic(Real a, Real a_value, Real a_deriv, Real b, Real b_value, Real b_deriv) const
        {
          // secant fallback
          if (!use_cubic || a_value == invalid || b_value == invalid)
            return secant(a,a_deriv,b,b_deriv);

          // quadratic interpolation on [0,b-a]
          const Real v_0 = a_deriv * std::pow(b - a,3);
          const Real v_1 = -2. * (3. * (b_value - a_value) * (b - a)
              - (b_deriv + 2. * a_deriv) * std::pow(b - a,2));
          const Real v_2 = 3. * ((b_deriv + a_deriv) * (b - a) - 2. * (b_value - a_value));

          // quadratic interpolation of derivative: d_2 x^2 + d_1 x + d_0
          const Real denom = std::pow(b - a,3);
          const Real d_0 = (v_0 + v_1 * a + v_2 * std::pow(a,2)) / denom;
          const Real d_1 = (-v_1 - 2. * v_2 * a) / denom;
          const Real d_2 = v_2 / denom;

          const Real discriminant = d_1*d_1 - 4.*d_0*d_2;
          Real cubic1, cubic2;
          if (discriminant < 0.)
          {
            if (std::abs(discriminant / (4.*d_0) / std::max(a_deriv,b_deriv)) < cubic_threshold)
              cubic1 = cubic2 = -2.*d_0/d_1;
            else
              return invalid;
          }
          else
          {
            cubic1 = -2.*d_0/(d_1 + std::sqrt(discriminant));
            cubic2 = -2.*d_0/(d_1 - std::sqrt(discriminant));
          }

          if (cubic1 > cubic2)
            std::swap(cubic1,cubic2);

          const Real sec = secant(a,a_deriv,b,b_deriv);

          if (verbosity >= 2)
            std::cout << " [ " << std::setw(13) << std::scientific << a
              << " , " << std::setw(13) << std::scientific << b
              << " ] , secant: " << std::setw(13) << std::scientific << sec
              << " , cubic: " << std::setw(13) << std::scientific << cubic1
              << " , " << std::setw(13) << std::scientific << cubic2 << std::endl;

          Real out;
          if (a_deriv < 0 && b_deriv < 0)
          {
            if (cubic2 < b)
              out = invalid;
            else if (cubic1 < b)
              out = cubic2;
            else
              out = cubic1;
          }
          else if (a_deriv > 0 && b_deriv > 0)
          {
            if (cubic1 > a)
              out = invalid;
            else if (cubic2 > a)
              out = cubic1;
            else
              out = cubic2;
          }
          else
          {
            const Real mid = 0.5 * (a + b);
            if (std::abs(cubic1 - mid) < std::abs(cubic2 - mid))
              out = cubic1;
            else
              out = cubic2;
          }

          return std::isfinite(out) ? out : invalid;
        }

        /**
         * @brief "secant^2" routine of CG_DESCENT line search, optionally with cubic interpolation
         *
         * Performs a secant interpolation step to shrink the search interval,
         * and then a secant extrapolation step from the part of the interval
         * that was discarded to further shrink the interval. In most cases,
         * this will result in both interval boundaries being updated.
         *
         * The secant interpolation / extrapolation steps are replaced with
         * cubic Hermite interpolation if use_cubic is true (default).
         *
         * @param         problem     problem definition, for function values and gradients
         * @param         direction   current search direction
         * @param         derivative  directional derivative at start of line
         * @param[in,out] point       current point, or updated point if accepted
         * @param[in,out] value       current value, or updated value if accepted
         * @param[in,out] gradient    current gradient, or updated gradient if accepted
         * @param[out]    deriv_alpha directional derivative at accepted step width
         * @param[in,out] a           left interval boundary
         * @param[in,out] a_value     value at left boundary
         * @param[in,out] a_deriv     derivative at left boundary
         * @param[in,out] b           right interval boundary
         * @param[in,out] b_value     value at right boundary
         * @param[in,out] b_deriv     derivative at right boundary
         *
         * @return true if acceptable step width was found, else false
         */
        bool doubleCubic(const ProblemBase<Real,Point>& problem, const Point& direction,
            Real derivative, Point& point, Real& value, Point& gradient, Real& deriv_alpha,
            Real& a, Real& a_value, Real& a_deriv, Real& b, Real& b_value, Real& b_deriv) const
        {
          const Real old_a = a;
          const Real old_a_value = a_value;
          const Real old_a_deriv = a_deriv;
          const Real old_b = b;
          const Real old_b_value = b_value;
          const Real old_b_deriv = b_deriv;

          Real c = cubic(a,a_value,a_deriv,b,b_value,b_deriv);

          if (updateInterval(problem,direction,derivative,point,value,gradient,
                deriv_alpha,a,a_value,a_deriv,b,b_value,b_deriv,c))
            return true;
          else if (c == b)
          {
            c = cubic(b,b_value,b_deriv,old_b,old_b_value,old_b_deriv);

            return updateInterval(problem,direction,derivative,point,value,gradient,
                deriv_alpha,a,a_value,a_deriv,b,b_value,b_deriv,c);
          }
          else if (c == a)
          {
            c = cubic(old_a,old_a_value,old_a_deriv,a,a_value,a_deriv);

            return updateInterval(problem,direction,derivative,point,value,gradient,
                deriv_alpha,a,a_value,a_deriv,b,b_value,b_deriv,c);
          }

          return false;
        }

        /**
         * @brief "update" routine of CG_DESCENT line search
         *
         * Updates the interval boundaries using information from the
         * secant / cubic interpolation steps, while guaranteeing that
         * the certain properties of the interval boundaries don't
         * change: basically, the right boundary needs to have positive
         * slope, the left boundary needs to have negative slope, and
         * the value at the left boundary needs to be smaller or equal
         * to the value from the previous step.
         *
         * If the result of the interpolation step fulfills one of these
         * properties, it replaces the corresponding boundary. If neither
         * of the two boundaries can be replaced, then a bisection
         * routine is carried out, starting with the interpolation
         * result and contracting the interval until a boundary is
         * successfully replaced.
         *
         * @param         problem     problem definition, for function values and gradients
         * @param         direction   current search direction
         * @param         derivative  directional derivative at start of line
         * @param[in,out] point       current point, or updated point if accepted
         * @param[in,out] value       current value, or updated value if accepted
         * @param[in,out] gradient    current gradient, or updated gradient if accepted
         * @param[out]    deriv_alpha directional derivative at accepted step width
         * @param[in,out] a           left interval boundary
         * @param[in,out] a_value     value at left boundary
         * @param[in,out] a_deriv     derivative at left boundary
         * @param[in,out] b           right interval boundary
         * @param[in,out] b_value     value at right boundary
         * @param[in,out] b_deriv     derivative at right boundary
         * @param         c           proposal for interval boundary update
         *
         * @return true if acceptable step width was found, else false
         */
        bool updateInterval(const ProblemBase<Real,Point>& problem, const Point& direction,
            Real derivative, Point& point, Real& value, Point& gradient, Real& deriv_alpha,
            Real& a, Real& a_value, Real& a_deriv, Real& b, Real& b_value, Real& b_deriv,
            const Real& c) const
        {
          if (c < a || c > b)
            return false;

          test_point = point;
          test_point.axpy(direction,c);

          if (combined_eval)
            std::tie(test_value, test_deriv) = problem.value_and_deriv(test_point,direction,test_grad);
          else
          {
            test_deriv = problem.dir_deriv(test_point,direction,test_grad);
            test_value = invalid;
          }

          if (verbosity >= 2)
            std::cout << "shrink (1) interval: [ " << std::setw(13) << std::scientific << a
              << " , " << std::setw(13) << std::scientific << b
              << " ] , width: " << std::setw(13) << std::scientific << b - a
              << " , probe: " << std::setw(13) << std::scientific << c
              << " , deriv: " << std::setw(13) << std::scientific << test_deriv << std::endl;

          if (checkCriterion(problem,direction,derivative,point,value,gradient,c,deriv_alpha))
          {
            a = c;
            b = c;
            return true;
          }

          // new point better than b, shrink from right
          if (test_deriv >= 0.)
          {
            b = c;
            b_value = test_value;
            b_deriv = test_deriv;
            return false;
          }
          else
          {
            if (test_value == invalid)
              test_value = problem.value(test_point,true);

            // new point better than a, shrink from left
            if (test_value <= value + epsilon * std::abs(value))
            {
              a = c;
              a_value = test_value;
              a_deriv = test_deriv;
              return false;
            }
            // better than neither, start bisection routine
            else
            {
              b = c;
              b_value = test_value;
              b_deriv = test_deriv;
              Real d = (1 - theta) * a + theta * b;
              if (verbosity >= 2)
                std::cout << "update (1) d = " << d << std::endl;

              test_point = point;
              test_point.axpy(direction,d);

              if (combined_eval)
                std::tie(test_value, test_deriv) = problem.value_and_deriv(test_point,direction,test_grad);
              else
              {
                test_deriv = problem.dir_deriv(test_point,direction,test_grad);
                test_value = invalid;
              }

              if (verbosity >= 2)
                std::cout << "shrink (2) interval: [ " << std::setw(13) << std::scientific << a
                  << " , " << std::setw(13) << std::scientific << b
                  << " ] , width: " << std::setw(13) << std::scientific << b - a
                  << " , probe: " << std::setw(13) << std::scientific << d
                  << " , deriv: " << std::setw(13) << std::scientific << test_deriv << std::endl;

              if (checkCriterion(problem,direction,derivative,point,value,gradient,d,deriv_alpha))
              {
                a = d;
                b = d;
                return true;
              }

              std::size_t count = 1;
              while (test_deriv < 0.)
              {
                if (test_value == invalid)
                  test_value = problem.value(test_point,true);

                if (test_value <= value + epsilon * std::abs(value))
                {
                  a = d;
                  a_value = test_value;
                  a_deriv = test_deriv;
                }
                else
                {
                  b = d;
                  b_value = test_value;
                  b_deriv = test_deriv;
                }

                d = (1 - theta) * a + theta * b;
                if (verbosity >= 2)
                  std::cout << "update (2) d = " << d << std::endl;
                test_point = point;
                test_point.axpy(direction,d);

                if (combined_eval)
                  std::tie(test_value, test_deriv) = problem.value_and_deriv(test_point,direction,test_grad);
                else
                {
                  test_deriv = problem.dir_deriv(test_point,direction,test_grad);
                  test_value = invalid;
                }

                if (verbosity >= 2)
                  std::cout << "shrink (3) interval: [ " << std::setw(13) << std::scientific << a
                    << " , " << std::setw(13) << std::scientific << b
                    << " ] , width: " << std::setw(13) << std::scientific << b - a
                    << " , probe: " << std::setw(13) << std::scientific << d
                    << " , deriv: " << std::setw(13) << std::scientific << test_deriv << std::endl;

                if (checkCriterion(problem,direction,derivative,point,value,gradient,d,deriv_alpha))
                {
                  a = d;
                  b = d;
                  return true;
                }

                count++;
                if (count > linesearch_tries)
                {
#if HAVE_DUNE_COMMON
                  DUNE_THROW(BracketingFailed,"failed to update linesearch interval");
#else
                  throw BracketingFailed();
#endif // HAVE_DUNE_COMMON
                }
              }

              b = d;
              b_value = test_value;
              b_deriv = test_deriv;
              return false;
            }
          }
        }

        /**
         * @brief Accept step width if line search criterion is fulfilled
         *
         * Checks the line search criterion for the suggested step width,
         * first the part depending on the directional derivative only, and
         * then the other parts. If the criterion is fulfilled, then the
         * current step width is accepted, and the old point with associated
         * function value, gradient, etc. is replaced by the new one.
         *
         * @param         problem     problem definition, for function values and gradients
         * @param         direction   current search direction
         * @param         derivative  directional derivative at start of line
         * @param[in,out] point       current point, or updated point if accepted
         * @param[in,out] value       current value, or updated value if accepted
         * @param[in,out] gradient    current gradient, or updated gradient if accepted
         * @param         alpha       suggested step width
         * @param[out]    deriv_alpha directional derivative at accepted step width
         *
         * @return whether line search criterion is fulfilled
         */
        bool checkCriterion(const ProblemBase<Real,Point>& problem,
            const Point& direction, Real derivative, Point& point, Real& value,
            Point& gradient, Real alpha, Real& deriv_alpha) const
        {
          if (this->criterion->check_derivative(value,derivative,alpha,test_deriv))
          {
            if (test_value == invalid)
              test_value = problem.value(test_point,true);

            if (this->criterion->check(value,derivative,alpha,test_value,test_deriv))
            {
              if (verbosity >= 1)
                std::cout << "Linesearch criterion satisfied, alpha: " << alpha << std::endl;

              problem.gradient(test_point,test_grad,&direction,test_deriv);

              std::swap(point,test_point);
              value = test_value;
              std::swap(gradient,test_grad);
              deriv_alpha = test_deriv;

              return true;
            }
          }

          return false;
        }

      };

  }
}

#endif // DUNE_NONLINOPT_LINESEARCH_HH
